rm(list=ls())

require(data.table)
# This scripts tries to get the table for analyzing the synthetic datasets
# d_cart = fread("./res_SIDES_and_CART.csv")
# # remove the columns if it contains "rsides"
# d_cart = d_cart[, !grepl("rsides", names(d_cart)), with=FALSE]
# # d_cart = fread("./res_cart_new.csv")
d_cart = fread("./res_cart_CameraReady.csv")
# rename the last four columns as                                      
#. gt_te_of_gt_subgroup = gt_te_of_gt_subgroup_,gt_te_of_learned_subgroup, estimated_te_of_learned_subgroup, estimated_te_of_gt_subgroup 
colnames(d_cart)[(ncol(d_cart)-3):ncol(d_cart)] = 
  c("gt_te_of_gt_subgroup", "estimated_te_of_gt_subgroup", 
    "gt_te_of_learned_subgroup", "estimated_te_of_learned_subgroup")
# get the summary table with each row represents the combination of simulator_name and n
#  collect the columns as follows: 
#. - mean jaccard plus/minus standard error
#  - mean tree_diff_te plus/minus standard error

diff_gt_te = abs(d_cart$gt_te_of_gt_subgroup - d_cart$gt_te_of_learned_subgroup)
diff_te_GtGt_EstimateLearned = abs(d_cart$gt_te_of_gt_subgroup - d_cart$estimated_te_of_learned_subgroup)
d_cart$diff_gt_te = diff_gt_te
d_cart$diff_te_GtGt_EstimateLearned = diff_te_GtGt_EstimateLearned
N_ = table(d_cart$n, d_cart$simulator_name)[1]
d_cart_summary = d_cart[, .(
  jaccard_mean = mean(tree_jaccard),
  jaccard_se = sd(tree_jaccard) / sqrt(N_),
  diff_te_mean = mean(tree_diff_te),
  diff_te_se = sd(tree_diff_te) / sqrt(N_),
  diff_te_sd = sd(tree_diff_te),
  jaccard_sd = sd(tree_jaccard), 
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]

######### RSIDES ########

d_rsides = fread("./res_SIDES_CameraReady.csv")
d_rsides$diff_gt_te = abs(d_rsides$gt_te_of_gt_subgroup - d_rsides$gt_te_of_learned_subgroup)
d_rsides$diff_te_GtGt_EstimateLearned = abs(d_rsides$gt_te_of_gt_subgroup - d_rsides$estimated_te_of_learned_subgroup)
# get the summary table with each row represents the combination of simulator_name and n
d_rsides_summary = d_rsides[, .(
  jaccard_mean = mean(rsides_jaccard),
  jaccard_se = sd(rsides_jaccard) / sqrt(N_),
  diff_te_mean = mean(rsides_diff_te),
  diff_te_se = sd(rsides_diff_te) / sqrt(N_),
  diff_te_sd = sd(rsides_diff_te),
  jaccard_sd = sd(rsides_jaccard),
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]

# # check 
# tmp1 = mean(d_cart_summary$jaccard_mean > d_rsides_summary$jaccard_mean)
# tmp2 = mean(d_cart_summary$diff_te_mean < d_rsides_summary$diff_te_mean)
# cat("Proportion of CART jaccard > RSIDES jaccard: ", tmp1, "\n",
#     "Proportion of CART tree_diff_te < RSIDES tree_diff_te: ", tmp2, "\n")


########### QUINT ##########
d_quint = fread("./res_quint_simulation_CameraReady.csv")
colnames(d_quint)[(ncol(d_quint)-3):ncol(d_quint)] = c("gt_te_of_gt_subgroup", "gt_te_of_learned_subgroup", 
                                                     "estimated_te_of_learned_subgroup", "estimated_te_of_gt_subgroup")

d_quint$diff_gt_te = abs(d_quint$gt_te_of_gt_subgroup - d_quint$gt_te_of_learned_subgroup)
d_quint$diff_te_GtGt_EstimateLearned = abs(d_quint$gt_te_of_gt_subgroup - d_quint$estimated_te_of_learned_subgroup)
d_quint_summary = d_quint[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard),
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_quint_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_quint_summary$diff_te_mean)
cat("Proportion of CART jaccard > QUINT jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < QUINT tree_diff_te: ", tmp2, "\n")

########### vt #########
# d_vt = fread("./res_vt_simulation.csv")
d_vt = fread("./res_vt_simulation_CR.csv")
colnames(d_vt)[(ncol(d_vt)-3):ncol(d_vt)] = c("gt_te_of_gt_subgroup", "gt_te_of_learned_subgroup", 
                                                   "estimated_te_of_learned_subgroup", "estimated_te_of_gt_subgroup")

d_vt$diff_gt_te = abs(d_vt$gt_te_of_gt_subgroup - d_vt$gt_te_of_learned_subgroup)
d_vt$diff_te_GtGt_EstimateLearned = abs(d_vt$gt_te_of_gt_subgroup - d_vt$estimated_te_of_learned_subgroup)
d_vt_summary = d_vt[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard),
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]

# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_vt_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_vt_summary$diff_te_mean)
cat("Proportion of CART jaccard > VT jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < VT tree_diff_te: ", tmp2, "\n")



######## Causal Tree ########
d_ct = fread("./competitors/causal_tree/res_CausalTree_simulation_CameraReady.csv")

d_ct$diff_gt_te = abs(d_ct$gt_te_of_gt_subgroup - d_ct$gt_te_of_learned_subgroup)
d_ct$diff_te_GtGt_EstimateLearned = abs(d_ct$gt_te_of_gt_subgroup - d_ct$estimated_te_of_learned_subgroup)


d_ct_summary = d_ct[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard),
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_ct_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_ct_summary$diff_te_mean)
cat("Proportion of CART jaccard > Causal Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Causal Tree tree_diff_te: ", tmp2, "\n")


####### Interaction Tree ######
d_it = fread("./competitors/interaction_tree/R code for IT_NOV2024/res_InteractionTree_simulation_CameraReady.csv")
d_it$diff_gt_te = abs(d_it$gt_te_of_gt_subgroup - d_it$gt_te_of_learned_subgroup)
d_it$diff_te_GtGt_EstimateLearned = abs(d_it$gt_te_of_gt_subgroup - d_it$estimated_te_of_learned_subgroup)

d_it_summary = d_it[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard),
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]

# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_it_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_it_summary$diff_te_mean)
cat("Proportion of CART jaccard > Interaction Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Interaction Tree tree_diff_te: ", tmp2, "\n")

####### distill tree
d_distill = fread("./competitors/distill_tree/res_distill_CameraReady.csv")
colnames(d_distill)[4] = "jaccard"
colnames(d_distill)[5] = "max_te"
colnames(d_distill)[7] = "diff_te"
d_distill$diff_gt_te = abs(d_distill$gt_te_of_gt_subgroup - d_distill$gt_te_of_learned_subgroup)
d_distill$diff_te_GtGt_EstimateLearned = abs(d_distill$gt_te_of_gt_subgroup - d_distill$estimated_te_of_learned_subgroup)

d_distill_summary = d_distill[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard), 
  diff_gt_te_mean = mean(diff_gt_te),
  diff_gt_te_se = sd(diff_gt_te) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_distill_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_distill_summary$diff_te_mean)
cat("Proportion of CART jaccard > Distill Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Distill Tree tree_diff_te: ", tmp2, "\n")

# CURLS
csv_files = list.files("./competitors/curls/res_server_cameraReady/myresults/synthetic/", full.names = TRUE, pattern = "*.csv")
data_info = basename(csv_files)
simulator_names = sapply(data_info, function(a){
  strsplit(a, split = "_n")[[1]][1]
})
simulator_names = unname(simulator_names)
require(dplyr)
n_values = sapply(data_info, function(a){
  unname(strsplit(a, split = "_n_")[[1]][2])
})
n_values = sapply(n_values, function(a){
  as.numeric(unname(strsplit(a, split = "_iter")[[1]][1]))
})
n_values = unname(n_values)
iter_values = sapply(data_info, function(a){
  unname(strsplit(a, split = "_iter_")[[1]][2])
})
iter_values = sapply(iter_values, function(a){
  as.numeric(unname(strsplit(a, split = "_curls")[[1]][1]))
})
iter_values = unname(iter_values)

d_curls = do.call(rbind, lapply(csv_files, read.csv))

# d_curls$diff_gt_estimated = abs(d_curls$diff_gt_estimated)
d_curls$simulator_name = simulator_names
d_curls$n = n_values
d_curls$iter = iter_values
d_curls = d_curls[d_curls$n != 500, ]

d_curls$diff_gt_estimated = abs(d_curls$gt_te_of_gt_subgroup - d_curls$gt_te_of_learned_subgroup)
d_curls$diff_te_GtGt_EstimateLearned = abs(d_curls$gt_te_of_gt_subgroup - d_curls$estimated_te_of_learned_subgroup)

d_curls = as.data.table(d_curls)
# order d_curls the rows by simulator_name and n, and the order of simulator_name follows d_cart
d_curls = d_curls[order(factor(simulator_name, levels = unique(d_cart$simulator_name)), n, iter)]

d_curls_summary = d_curls[, .(
  jaccard_mean = mean(jaccard_similarity),
  jaccard_se = sd(jaccard_similarity) / sqrt(N_),
  diff_te_mean = mean(diff_gt_estimated),
  diff_te_se = sd(diff_gt_estimated) / sqrt(N_),
  diff_te_sd = sd(diff_gt_estimated),
  jaccard_sd = sd(jaccard_similarity),
  diff_gt_te_mean = mean(diff_gt_estimated),
  diff_gt_te_se = sd(diff_gt_estimated) / sqrt(N_),
  diff_te_GtGt_EstimateLearned_mean = mean(diff_te_GtGt_EstimateLearned),
  diff_te_GtGt_EstimateLearned_se = sd(diff_te_GtGt_EstimateLearned) / sqrt(N_)
), by = .(simulator_name, n)]
# remove the rows with n == 500
# re-order the rows by simulator_name and n, and the order of simulator_name follows d_cart
d_curls_summary = d_curls_summary[order(factor(simulator_name, levels = unique(d_cart_summary$simulator_name)), n,)]

# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_curls_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_curls_summary$diff_te_mean)
cat("Proportion of CART jaccard > CURLS jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < CURLS tree_diff_te: ", tmp2, "\n")


# combine all the summaries, and add a column for "algorithm"
d_combined_summary = rbind(
  d_cart_summary[, algorithm := "Ours"],
  d_rsides_summary[, algorithm := "SIDES"],
  d_quint_summary[, algorithm := "QUINT"],
  d_vt_summary[, algorithm := "VT"],
  d_ct_summary[, algorithm := "Causal Tree"],
  d_it_summary[, algorithm := "Interaction Tree"],
  d_curls_summary[, algorithm := "CURLS"],
  d_distill_summary[, algorithm := "Distill Tree"]
)

# replace the 3 simulator_name with "simulation1, simulation2, ..."
d_combined_summary$simulator_name = factor(d_combined_summary$simulator_name,
                                           levels = unique(d_combined_summary$simulator_name),
                                           labels = paste0("Simulation ", seq_along(unique(d_combined_summary$simulator_name))))
d_combined_summary$algorithm <- factor(
  d_combined_summary$algorithm,
  levels = c("CURLS", "SIDES", "VT", "Distill Tree", "QUINT", "Causal Tree", "Interaction Tree", "Ours")
)

# make figures
library(ggplot2)
library(patchwork)
library(data.table)


###########
# for diff_gt_estimated
##########
setDT(d_combined_summary)
sim_names <- unique(d_combined_summary$simulator_name)
plots <- vector("list", length(sim_names))

require(RColorBrewer)
base_cols <- brewer.pal(7, "Dark2")
names(base_cols) <- c("CURLS","VT","QUINT","Interaction Tree",
                      "SIDES","Distill Tree","Causal Tree")

# 2. append “Ours” as black
all_cols <- c(base_cols, Ours = "black")

for (i in seq_along(sim_names)) {
  sim_name <- sim_names[i]
  
  p <- ggplot(d_combined_summary[simulator_name == sim_name], 
              aes(x = n, y = diff_gt_te_mean, color = algorithm)) +
    geom_point() +
    geom_line() +
    scale_size_identity() + 
    scale_linewidth_identity() +
    geom_errorbar(aes(ymin = diff_gt_te_mean - diff_gt_te_se, 
                      ymax = diff_gt_te_mean + diff_gt_te_se), width = 0.2) +
    labs(title = sim_name) + 
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_manual(values = all_cols)
  
  # ylim(0, max(d_combined_summary[simulator_name == sim_name]$diff_te_mean + d_combined_summary[simulator_name == sim_name]$diff_te_se) + 0.05) + 
  
  # Add x-axis label ONLY to the last plot
  if (i == length(sim_names) - 1) {
    p <- p + labs(x = "Sample Size (n)") +
      theme(axis.title.x = element_text(), 
            legend.position = "bottom", 
            legend.title = element_blank())  # force it to show
  }
  
  plots[[i]] <- p
}

# Combine plots
final_plot <- wrap_plots(plots, ncol = 3) &
  theme(axis.title.y = element_text(angle = 90)) &
  labs(y = "Mean Absolute Error \n of Ground-Truth Treatment Effect \n of the Learned vs. Ground-Truth Subgroup")
final_plot[[2]] = final_plot[[2]] + 
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot[[length(plots)]] = final_plot[[length(plots)]] +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot
ggsave("./CR_figures/Camera_diff_GtGt_GtLearned.png", plot = final_plot, width = 6, height = 5, dpi = 300)



#######################
##### For Diff TE #####
#######################

setDT(d_combined_summary)
sim_names <- unique(d_combined_summary$simulator_name)
plots <- vector("list", length(sim_names))

require(RColorBrewer)
base_cols <- brewer.pal(7, "Dark2")
names(base_cols) <- c("CURLS","VT","QUINT","Interaction Tree",
                      "SIDES","Distill Tree","Causal Tree")

# 2. append “Ours” as black
all_cols <- c(base_cols, Ours = "black")

for (i in seq_along(sim_names)) {
  sim_name <- sim_names[i]
  
  p <- ggplot(d_combined_summary[simulator_name == sim_name], 
              aes(x = n, y = diff_te_GtGt_EstimateLearned_mean, color = algorithm)) +
    geom_point() +
    geom_line() +
    scale_size_identity() + 
    scale_linewidth_identity() +
    geom_errorbar(aes(ymin = diff_te_GtGt_EstimateLearned_mean - diff_te_GtGt_EstimateLearned_se, 
                      ymax = diff_te_GtGt_EstimateLearned_mean + diff_te_GtGt_EstimateLearned_se), width = 0.2) +
    labs(title = sim_name) + 
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_manual(values = all_cols)
  
  # ylim(0, max(d_combined_summary[simulator_name == sim_name]$diff_te_mean + d_combined_summary[simulator_name == sim_name]$diff_te_se) + 0.05) + 
  
  # Add x-axis label ONLY to the last plot
  if (i == length(sim_names) - 1) {
    p <- p + labs(x = "Sample Size (n)") +
      theme(axis.title.x = element_text(), 
            legend.position = "bottom", 
            legend.title = element_blank())  # force it to show
  }
  
  plots[[i]] <- p
}

# Combine plots
final_plot <- wrap_plots(plots, ncol = 3) &
  theme(axis.title.y = element_text(angle = 90)) &
  labs(y = "Mean Absolute Error of \n Estimated Treatment Effect \n of The Learned vs. Ground-Truth Subgroup")
final_plot[[2]] = final_plot[[2]] + 
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot[[length(plots)]] = final_plot[[length(plots)]] +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot

ggsave("./CR_figures/Camera_diff_te_GtGt_EstimateLearned_combined.png", plot = final_plot, width = 6, height = 5, dpi = 300)

#######################
##### For Jaccard #####
#######################

setDT(d_combined_summary)
sim_names <- unique(d_combined_summary$simulator_name)
plots <- vector("list", length(sim_names))

for (i in seq_along(sim_names)) {
  sim_name <- sim_names[i]
  
  p <- ggplot(d_combined_summary[simulator_name == sim_name], 
              aes(x = n, y = jaccard_mean, color = algorithm)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = jaccard_mean - jaccard_se, ymax = jaccard_mean + jaccard_se), width = 0.2) +
    labs(title = sim_name) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_manual(values = all_cols)
  
  # ylim(0, max(d_combined_summary[simulator_name == sim_name]$jaccard_mean + d_combined_summary[simulator_name == sim_name]$jaccard_se) + 0.05) + 
  
  # Add x-axis label ONLY to the last plot
  if (i == length(sim_names) - 1) {
    p <- p + labs(x = "Sample Size (n)") +
      theme(axis.title.x = element_text(), 
            legend.position = "bottom", 
            legend.title = element_blank())  # force it to show
  }
  
  plots[[i]] <- p
}

# Combine plots
final_plot <- wrap_plots(plots, ncol = 3) &
  theme(axis.title.y = element_text(angle = 90, vjust = 1)) &
  labs(y = "Jaccard Similarity with\n The Ground Truth Subgroup")
final_plot[[2]] = final_plot[[2]] + 
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot[[length(plots)]] = final_plot[[length(plots)]] +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "")
final_plot


ggsave("./CR_figures/Camera_jaccard_combined.png", plot = final_plot, width = 6, height = 5, dpi = 300)

# check the difference between the d_combined_summary$diff_te_mean for CART and Distill
d_combined_summary[(algorithm == "Ours" | algorithm == "Distill Tree") & n == 5000]
     