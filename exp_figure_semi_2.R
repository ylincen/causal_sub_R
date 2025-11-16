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

d_combine2 = d_combine[d_combine$simulator == "ACIC2016",]

d_combine2$diff_with_gt_ratio = abs(d_combine2$max_subgroup_diff_average) / 
  abs(d_combine2$max_subgroup_gt_treatment_effct)
d_combine2$estimated_as_ratio = d_combine2$max_subgroup_treatment_effect / 
  d_combine2$max_subgroup_gt_treatment_effct

d_combine2$dataset_index = sapply(d_combine2$dataset, function(a){
  as.numeric(sub("ACIC2016_(\\d+)\\.csv", "\\1", a))
})
# re-order the rows by the index
# d_combine2 = d_combine2[order(d_combine2$dataset_index), ]


# boxplot for accuracy

d_combine2$algorithm = factor(d_combine2$algorithm, 
                              levels = c("Ours", "SIDES", "QUINT", 
                                         "CausalTree", "InteractionTree", 
                                         "DistillTree", "CURLS"))

# get the average max_subgroup_diff_average by algorithm, using the aggregation of data.table
d_combine2[,
  mean(max_subgroup_diff_average, na.rm = TRUE), 
  by = algorithm]


ggplot(d_combine2, aes(x = algorithm, 
                         y = max_subgroup_diff_average, 
                         fill = algorithm)) + 
  geom_boxplot() +
  ylab("Difference between the Estimated and Ground-truth Treatment Effect for the Subgroups Learned From Data") +
  xlab("Algorithm") +
  theme_bw() +
  theme(legend.position = "none")


d_combine2$dataset_label <- sub("ACIC2016_(\\d+)\\.csv", "\\1", d_combine2$dataset)
# ordered_ds <- d_combine2[order(dataset_index), unique(dataset)]
# # overwrite dataset as a factor
# d_combine2[, dataset := factor(dataset, levels = ordered_ds)]

for(i in 1:4){
  rows = ((i - 1) * 20 + 1):(i * 20)
  datasets_selected = paste0("ACIC2016_", rows, ".csv")
  d_combine2_part = d_combine2[d_combine2$dataset %in% datasets_selected, ]
  
  
  p = ggplot(d_combine2_part, aes(x=dataset, y=max_subgroup_gt_treatment_effct, fill=algorithm,
                         color=algorithm)) + 
    geom_point(aes(size = ifelse(algorithm == "Ours", 3, 2), 
                   alpha = ifelse(algorithm=="Ours", 1, 0.6),
                   shape = ifelse(algorithm=="Ours", "square", "round"))) +
    scale_size_identity() + 
    scale_alpha_identity() + 
    scale_x_discrete(name = "Semi-synthetic Datasets",
                     labels = setNames(d_combine2$dataset_label, d_combine2$dataset)) + 
    ylab("Ground Truth Treatment effect of Learned Subgroups") + 
    xlab("Semi-synthetic Dataset Index") + 
    guides(color = guide_legend(title = "algorithm"),
           linetype = "none", size = "none", shape = "none")
  # ggsave(paste0("./figures/semi_figures/gt_subgroups_part_", i, ".png"), p, 
  #        width = 6, height = 4, device = "png", dpi = 600)
}

p_list <- vector("list", 4)

for (i in 1:4) {
  rows <- ((i - 1) * 20 + 1):(i * 20)
  datasets_selected <- paste0("ACIC2016_", rows, ".csv")
  d_part <- subset(d_combine2, dataset %in% datasets_selected)
  sel  <- paste0("ACIC2016_", rows, ".csv")
  d_part[, dataset := factor(dataset, levels = sel)]
  
  
  p <- ggplot(d_part,
              aes(x = dataset,
                  y = max_subgroup_gt_treatment_effct,
                  fill = algorithm, color = algorithm)) +
    geom_point(aes(size  = ifelse(algorithm == "Ours", 3, 2),
                   alpha = ifelse(algorithm == "Ours", 1, 0.6))) +    # 15=square,16=circle
    scale_size_identity() +
    scale_alpha_identity() +
    scale_shape_identity() +
    scale_x_discrete(name = "Semi-synthetic Datasets",
                     labels = setNames(d_combine2$dataset_label, d_combine2$dataset)) +
    ylab("Average Treatment effect \n of Learned Subgroups \n obtained from the ground-\ntruth counterfactual outcomes") +
    xlab("Semi-synthetic Dataset Index") +
    guides(color = guide_legend(title = "algorithm"),
           linetype = "none", size = "none", shape = "none",alpha="none")
  # scale_x_discrete(name = "Semi-synthetic Datasets",
  #                  labels = setNames(d_combine2$dataset_label, d_combine2$dataset)) +
  p_list[[i]] <- p
}

# 2x2 layout, single (collected) legend at bottom, auto panel tags (A–D)
combined <- wrap_plots(plotlist = p_list, ncol = 1) +
  plot_layout(guides = "collect")&
  theme(legend.position = "bottom")

combined

ggsave("./figures/semi_figures/gt_subgroups_all.png",
       combined, width = 6, height = 10, dpi = 600)

# 
#   
# 
# # ggplot(d_combine2, aes(x=algorithm, y=max_subgroup_gt_treatment_effct, fill=algorithm,
# #                                                             color=algorithm)) + 
# #   geom_boxplot()
# 
# d_combine2$dataset = factor(d_combine2$dataset, 
#                             levels = paste0("ACIC2016_", 1:length(unique(d_combine2$dataset)), ".csv"))
# 
# d_combine2$lower = ifelse(d_combine2$max_subgroup_treatment_effect > 
#                     d_combine2$max_subgroup_gt_treatment_effct, 
#                     d_combine2$max_subgroup_gt_treatment_effct, 
#                     d_combine2$max_subgroup_treatment_effect)
# d_combine2$upper = ifelse(d_combine2$max_subgroup_treatment_effect <
#                     d_combine2$max_subgroup_gt_treatment_effct, 
#                     d_combine2$max_subgroup_gt_treatment_effct, 
#                     d_combine2$max_subgroup_treatment_effect)
# 
# # levels(d_combine2$dataset) = levels(d_combine2$dataset)[order(variance_te)]
# d_combine2$dataset_label <- sub("ACIC2016_(\\d+)\\.csv", "\\1", d_combine2$dataset)
# 
# for(alg in c("SIDES", "QUINT", "CausalTree", "InteractionTree", "DistillTree", "CURLS")){
#   p = ggplot(d_combine2[d_combine2$algorithm %in% c("Ours", alg)], 
#              aes(x = dataset, 
#                  y = max_subgroup_gt_treatment_effct, 
#                  color = algorithm,
#                  group = algorithm)) + 
#     geom_point(aes(size = ifelse(algorithm == "Ours", 3, 2))) +
#     scale_size_identity() +
#     scale_x_discrete(name = "Semi-synthetic Datasets",
#                      labels = setNames(d_combine2$dataset_label, d_combine2$dataset)) + 
#     ylab("max_subgroup_gt_treatment_effct") +
#     theme_bw()
#   dir = "./figures/semi_figures/"
#   if(!dir.exists(dir)){
#     dir.create(dir, recursive = TRUE)
#   }
#   # ggsave(paste0("./figures/semi_figures/semi_synthetic_ours_vs_", alg, ".png"), p, width = 8, height = 4,
#   #        device = "png", dpi = 600)
# }
# 
# d_combine2_wide = dcast(d_combine2, dataset ~ algorithm, value.var = "max_subgroup_gt_treatment_effct")
# for(baselines in c("SIDES", "QUINT", "CausalTree", "InteractionTree", "DistillTree", "CURLS")){
#   d_combine2_wide[[paste0("subgroup_gt_ours_vs_", baselines)]] = d_combine2_wide$Ours - d_combine2_wide[[baselines]]
#   cat("Percentage of datasets where Ours is better than ", baselines, ": ",
#       mean(d_combine2_wide[[paste0("subgroup_gt_ours_vs_", baselines)]] > 0, na.rm = TRUE), "\n")
# }
# 
# d_combine2_wide_vs = d_combine2_wide[, c("dataset",
#                                         "subgroup_gt_ours_vs_SIDES", 
#                                         "subgroup_gt_ours_vs_QUINT",
#                                         "subgroup_gt_ours_vs_CausalTree",
#                                         "subgroup_gt_ours_vs_InteractionTree",
#                                         "subgroup_gt_ours_vs_DistillTree",
#                                         "subgroup_gt_ours_vs_CURLS")]
# # to long
# d_combine2_vs = melt(d_combine2_wide_vs, id.vars = "dataset", 
#                                  variable.name = "comparison", 
#                                  value.name = "subgroup_gt_diff_ours_vs_baseline")
# ggplot(d_combine2_vs, aes(x = comparison, y= subgroup_gt_diff_ours_vs_baseline, 
#                           fill = comparison)) + 
#   geom_boxplot()
# 
# 
# d_combine2_wide_acc = dcast(d_combine2, dataset ~ algorithm, value.var = "max_subgroup_diff_average")
# for(baselines in c("SIDES", "QUINT", "CausalTree", "InteractionTree", "DistillTree", "CURLS")){
#   d_combine2_wide_acc[[paste0("subgroup_diff_ours_vs_", baselines)]] = d_combine2_wide_acc$Ours <= d_combine2_wide_acc[[baselines]]
#   cat("Percentage of datasets where Ours is better than ", baselines, ": ",
#       mean(d_combine2_wide_acc[[paste0("subgroup_diff_ours_vs_", baselines)]], na.rm = TRUE), "\n")
# }
# 
# 
# 
# 
# ggplot(d_combine2[d_combine2$algorithm %in% c("Ours", "CausalTree")], aes(x = dataset,
#                         y = max_subgroup_gt_treatment_effct,
#                         color = algorithm,
#                         group = algorithm)) +
#   geom_point(aes(size = ifelse(algorithm == "Ours", 3, 2))) +
#   geom_errorbar(aes(ymin = lower,
#                         ymax = upper),
#                 width = 1) +
#   scale_size_identity() +
#   scale_x_discrete(name = "Semi-synthetic Datasets",
#                    labels = setNames(d_combine2$dataset_label, d_combine2$dataset)) +
#   ylab("Max Subgroup Difference (Average)") +
#   theme_bw()
# 
# 
# 
# tmp1 = d_combine2[, c("dataset", "max_subgroup_gt_treatment_effct", "algorithm")]
# tmp1$type = "ground-truth"
# tmp2 = d_combine2[, c("dataset", "max_subgroup_treatment_effect", "algorithm")]
# tmp2$type = "estimated"
# colnames(tmp1)[2] = "treatment_effect"
# colnames(tmp2)[2] = "treatment_effect"
# d_combine3 = rbind(tmp1,
#                    tmp2)
# 
# ggplot(d_combine3, aes(x = dataset,
#                         y = treatment_effect,
#                         color = algorithm,
#                         shape = type)) +
#   geom_point(aes(size = ifelse(algorithm == "Ours", 3, 2))) +
#   scale_size_identity() +
#   scale_x_discrete(name = "Semi-synthetic Datasets") +
#   ylab("Treatment Effect") +
#   theme_bw() + ylim(c(0,20))
# 
# #
# #
# # d_combine2$diff_percetange = abs(d_combine2$max_subgroup_diff_average - d_combine2$max_subgroup_gt_treatment_effct) /
# #   abs(d_combine2$max_subgroup_gt_treatment_effct)
# # # ggplot(d_combine2, aes(x = dataset, y = max_subgroup_gt_treatment_effct, color=algorithm,
# # #                        group = algorithm)) +
# # #   geom_point()
# # ggplot(d_combine2, aes(x = as.numeric(factor(dataset)),
# #                        y = max_subgroup_gt_treatment_effct,
# #                        color = algorithm,
# #                        group = algorithm)) +
# #   geom_point(aes(size = ifelse(algorithm == "Ours", 4, 2))) +
# #   scale_size_identity() +
# #   scale_x_continuous(breaks = 1:length(unique(d_combine2$dataset)),
# #                      labels = unique(d_combine2$dataset),
# #                      name = "Semi-synthetic Datasets") +
# #   ylab("Ground Truth Treatment Effect") +
# #   theme_bw()
# #
# # ggplot(d_combine2, aes(x = as.numeric(factor(dataset)),
# #                        y = abs(max_subgroup_gt_treatment_effct - max_subgroup_treatment_effect) / abs(max_subgroup_gt_treatment_effct),
# #                        color = algorithm,
# #                        group = algorithm)) +
# #   geom_point(aes(size = ifelse(algorithm == "Ours", 4, 2))) +
# #   scale_size_identity() +
# #   scale_x_continuous(breaks = 1:length(unique(d_combine2$dataset)),
# #                      labels = unique(d_combine2$dataset),
# #                      name = "Semi-synthetic Datasets") +
# #   ylab("Ground Truth Treatment Effect") +
# #   ylim(c(0,1))
# #   theme_bw()
# #
# #
# # d_combine2_wide_acc = dcast(d_combine2, dataset ~ algorithm, value.var = "max_subgroup_diff_average")
# # apply(d_combine2_wide_acc[, 2:ncol(d_combine2_wide_acc)], 1, rank, na.last="keep") %>%
# #   apply(1, mean, na.rm=T)
# #
# # d_combine2_wide_acc_percentage = dcast(d_combine2, dataset ~ algorithm, value.var = "diff_percetange")
# # apply(d_combine2_wide_acc_percentage[, 2:ncol(d_combine2_wide_acc_percentage)], 1, rank, na.last="keep") %>%
# #   apply(1, mean, na.rm=T)
# #
# #
# # d_combine2_wide = dcast(d_combine2, dataset ~ algorithm, value.var = "max_subgroup_gt_treatment_effct")
# #
# # get_rank_percentage = function(d_combine2_wide, alg_name){
# #   res1 = mean(d_combine2_wide[[alg_name]] >= apply(d_combine2_wide[,2:ncol(d_combine2_wide)], 1, max, na.rm=T), na.rm=T)
# #   res2 = mean(d_combine2_wide[[alg_name]] >= apply(d_combine2_wide[,2:ncol(d_combine2_wide)], 1, sort, decreasing=T) %>%
# #          sapply(function(a) a[2]), na.rm=T)
# #   # print(c(res1, res2))
# #   return(c(res1, res2))
# # }
# #
# # for(alg_name in c("Ours", "SIDES", "QUINT", "CausalTree", "InteractionTree", "DistillTree", "CURLS")){
# #   cat("Algorithm: ", alg_name, "\n")
# #   res = get_rank_percentage(d_combine2_wide, alg_name)
# #   cat("Percentage of datasets where the algorithm is the best: ", res[1], "\n")
# #   cat("Percentage of datasets where the algorithm is the second best: ", res[2], "\n")
# # }
# #
# #
# #
# #
# #
