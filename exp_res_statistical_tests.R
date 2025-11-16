rm(list=ls())

require(data.table)
# This scripts tries to get the table for analyzing the synthetic datasets
d_cart = fread("./res_SIDES_and_CART.csv")
# remove the columns if it contains "rsides"
d_cart = d_cart[, !grepl("rsides", names(d_cart)), with=FALSE]
# d_cart = fread("./res_cart_new.csv")

# get the summary table with each row represents the combination of simulator_name and n
#  collect the columns as follows: 
#. - mean jaccard plus/minus standard error
#  - mean tree_diff_te plus/minus standard error

N_ = table(d_cart$n, d_cart$simulator_name)[1]
d_cart_summary = d_cart[, .(
  jaccard_mean = mean(tree_jaccard),
  jaccard_se = sd(tree_jaccard) / sqrt(N_),
  diff_te_mean = mean(tree_diff_te),
  diff_te_se = sd(tree_diff_te) / sqrt(N_),
  diff_te_sd = sd(tree_diff_te),
  jaccard_sd = sd(tree_jaccard)
), by = .(simulator_name, n)]

######### RSIDES ########

d_rsides = fread("./res_SIDES_and_CART.csv")
d_rsides = d_rsides[, !grepl("cart", names(d_rsides)), with=FALSE]
# get the summary table with each row represents the combination of simulator_name and n
d_rsides_summary = d_rsides[, .(
  jaccard_mean = mean(rsides_jaccard),
  jaccard_se = sd(rsides_jaccard) / sqrt(N_),
  diff_te_mean = mean(rsides_diff_te),
  diff_te_se = sd(rsides_diff_te) / sqrt(N_),
  diff_te_sd = sd(rsides_diff_te),
  jaccard_sd = sd(rsides_jaccard)
), by = .(simulator_name, n)]

# check 
tmp1 = mean(d_cart_summary$jaccard_mean > d_rsides_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_rsides_summary$diff_te_mean)
cat("Proportion of CART jaccard > RSIDES jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < RSIDES tree_diff_te: ", tmp2, "\n")


########### QUINT ##########
d_quint = fread("./res_quint_simulation.csv")
d_quint_summary = d_quint[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_quint_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_quint_summary$diff_te_mean)
cat("Proportion of CART jaccard > QUINT jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < QUINT tree_diff_te: ", tmp2, "\n")

########### vt #########
d_vt = fread("./res_vt_simulation.csv")
d_vt_summary = d_vt[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard)
), by = .(simulator_name, n)]

# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_vt_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_vt_summary$diff_te_mean)
cat("Proportion of CART jaccard > VT jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < VT tree_diff_te: ", tmp2, "\n")



######## Causal Tree ########
d_ct = fread("./competitors/causal_tree/res_CausalTree_simulation.csv")
d_ct_summary = d_ct[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_ct_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_ct_summary$diff_te_mean)
cat("Proportion of CART jaccard > Causal Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Causal Tree tree_diff_te: ", tmp2, "\n")


####### Interaction Tree ######
d_it = fread("./competitors/interaction_tree/R code for IT_NOV2024/res_InteractionTree_simulation.csv")
d_it_summary = d_it[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard)
), by = .(simulator_name, n)]

# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_it_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_it_summary$diff_te_mean)
cat("Proportion of CART jaccard > Interaction Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Interaction Tree tree_diff_te: ", tmp2, "\n")

####### distill tree
d_distill = fread("./competitors/distill_tree/res_distill.csv")
colnames(d_distill)[4] = "jaccard"
colnames(d_distill)[5] = "max_te"
colnames(d_distill)[7] = "diff_te"

d_distill_summary = d_distill[, .(
  jaccard_mean = mean(jaccard),
  jaccard_se = sd(jaccard) / sqrt(N_),
  diff_te_mean = mean(diff_te),
  diff_te_se = sd(diff_te) / sqrt(N_),
  diff_te_sd = sd(diff_te),
  jaccard_sd = sd(jaccard)
), by = .(simulator_name, n)]
# check
tmp1 = mean(d_cart_summary$jaccard_mean > d_distill_summary$jaccard_mean)
tmp2 = mean(d_cart_summary$diff_te_mean < d_distill_summary$diff_te_mean)
cat("Proportion of CART jaccard > Distill Tree jaccard: ", tmp1, "\n",
    "Proportion of CART tree_diff_te < Distill Tree tree_diff_te: ", tmp2, "\n")

# CURLS
csv_files = list.files("./competitors/curls/res_server/synthetic/", full.names = TRUE, pattern = "*.csv")
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

d_curls$diff_gt_estimated = abs(d_curls$diff_gt_estimated)
d_curls$simulator_name = simulator_names
d_curls$n = n_values
d_curls$iter = iter_values
d_curls = d_curls[d_curls$n != 500, ]
d_curls = as.data.table(d_curls)
# order d_curls the rows by simulator_name and n, and the order of simulator_name follows d_cart
d_curls = d_curls[order(factor(simulator_name, levels = unique(d_cart$simulator_name)), n, iter)]

d_curls_summary = d_curls[, .(
  jaccard_mean = mean(jaccard_similarity),
  jaccard_se = sd(jaccard_similarity) / sqrt(N_),
  diff_te_mean = mean(diff_gt_estimated),
  diff_te_se = sd(diff_gt_estimated) / sqrt(N_),
  diff_te_sd = sd(diff_gt_estimated),
  jaccard_sd = sd(jaccard_similarity)
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

# long to wide
d_combined_summary_wide_jaccard_mean = 
  d_combined_summary[, .(simulator_name, n, algorithm, jaccard_mean)] %>%
  dcast(simulator_name + n ~ algorithm, value.var = "jaccard_mean")

# pairwise.wilcox.test with holm correction
library(data.table)

# long table: one row per (simulator_name, n, algorithm)
# columns: simulator_name, n, algorithm, jaccard_mean (or your metric)
dt <- copy(d_combined_summary)
score_col <- "jaccard_mean"   # change to your metric if needed
baseline  <- "Ours"
alt       <- "greater"        # "less" or "two.sided" if that’s your hypothesis

# Wide by dataset so we can pair by (simulator_name, n)
wide <- dcast(dt, simulator_name + n ~ algorithm, value.var = score_col)

# Keep datasets where baseline is present
stopifnot(baseline %in% names(wide))
methods <- setdiff(names(wide), c("simulator_name", "n", baseline))

res <- rbindlist(lapply(methods, function(m) {
  base <- wide[[baseline]]
  other <- wide[[m]]
  
  # matched pairs across datasets
  cc <- complete.cases(base, other)
  base <- base[cc]; other <- other[cc]
  
  # Wilcoxon ignores zero-diff pairs; handle all-zero diffs gracefully
  nz <- (base - other) != 0
  if (!any(nz)) {
    return(data.table(method = m, n_pairs = length(base), W = NA_real_, p = 1))
  }
  
  w <- wilcox.test(base[nz], other[nz], paired = TRUE, alternative = alt)
  
  data.table(method = m,
             n_pairs = sum(nz),
             W = unname(w$statistic),
             p = w$p.value)
}), use.names = TRUE, fill = TRUE)

# Adjust across all Ours-vs-X comparisons
res[, p_holm := p.adjust(p, method = "holm")][order(p_holm)]
res


# for diff_te_mean
dt <- copy(d_combined_summary)
score_col <- "diff_te_mean"   # change to your metric if needed
baseline  <- "Ours"
alt       <- "less"        # "less" or "two.sided" if that’s your hypothesis

# Wide by dataset so we can pair by (simulator_name, n)
wide <- dcast(dt, simulator_name + n ~ algorithm, value.var = score_col)

# Keep datasets where baseline is present
stopifnot(baseline %in% names(wide))
methods <- setdiff(names(wide), c("simulator_name", "n", baseline))

res <- rbindlist(lapply(methods, function(m) {
  base <- wide[[baseline]]
  other <- wide[[m]]
  
  # matched pairs across datasets
  cc <- complete.cases(base, other)
  base <- base[cc]; other <- other[cc]
  
  # Wilcoxon ignores zero-diff pairs; handle all-zero diffs gracefully
  nz <- (base - other) != 0
  if (!any(nz)) {
    return(data.table(method = m, n_pairs = length(base), W = NA_real_, p = 1))
  }
  
  w <- wilcox.test(base[nz], other[nz], paired = TRUE, alternative = alt)
  
  data.table(method = m,
             n_pairs = sum(nz),
             W = unname(w$statistic),
             p = w$p.value)
}), use.names = TRUE, fill = TRUE)

# Adjust across all Ours-vs-X comparisons
res[, p_holm := p.adjust(p, method = "holm")][order(p_holm)]
res

res$p <- signif(res$p, 4)
res$p_holm <- signif(res$p_holm, 4)
res$n_pairs = NULL
res$W = NULL
res
