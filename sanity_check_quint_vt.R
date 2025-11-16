rm(list = ls())

library(aVirtualTwins)
library(quint)
library(data.table)

source("./util_synthetic.R")

## ----------------------------------------------------------
## Setup
## ----------------------------------------------------------

simulator_names <- c(
  "simulate1",
  "simulate_imbalance_treatment",
  "simulate_long_rule"
)

simulator_names <- c(
  "simulate1"
)

n <- 1000
iters <- 1:1

## Helper: QUINT rule extraction (as in your code)
get_rules <- function(df) {
  kids <- strsplit(df$childnodes, ",")
  names(kids) <- df$parentnode
  
  rules <- list()
  walk_tree <- function(node, path = character()) {
    row <- match(node, df$parentnode)
    
    if (is.na(row)) {
      rules[[length(rules) + 1]] <<- paste(path, collapse = " & ")
      return(invisible())
    }
    
    v  <- df$splittingvar[row]
    sp <- df$truesplitpoint[row]
    ch <- as.integer(kids[[row]])
    
    walk_tree(ch[1], c(path, sprintf("%s < %s", v, sp)))
    walk_tree(ch[2], c(path, sprintf("%s >= %s", v, sp)))
  }
  
  walk_tree(1)
  unlist(rules, use.names = FALSE)
}

## Small helper: build hash of subgroup membership
subgroup_hash_fun <- function(bool_vec) {
  idx <- which(bool_vec)
  if (length(idx) == 0L) return("EMPTY")
  paste(idx, collapse = "-")
}

## ----------------------------------------------------------
## Main loop: run VT and QUINT, log best rules & TE
## ----------------------------------------------------------

rows_list <- list()
counter <- 0L

for (simulator_name in simulator_names) {
  for (iter_ in iters) {
    counter <- counter + 1L
    set.seed(counter)
    
    cat("Simulator:", simulator_name, "n =", n, "iter =", iter_, "\n")
    
    file_path <- paste0("./new_simulation/", simulator_name,
                        "/n_", n, "_iter_", iter_, ".csv")
    d <- read.csv(file_path)
    
    # train/test split
    train_indices <- sample(1:n, n * 0.5, replace = FALSE)
    d_train <- d[train_indices, ]
    d_test  <- d[-train_indices, ]
    
    # Ground truth quantities (shared by both methods)
    gt_bool <- get_gt_bool(dd = d_test, simulator_name = simulator_name)
    gt_te <- mean(d_test$Y[gt_bool & (d_test$T == 1)]) -
      mean(d_test$Y[gt_bool & (d_test$T == 0)])
    gt_te_theoretical_per_sample <- gt_te_per_sample(d_test, simulator_name)
    gt_te_of_gt_subgroup <- mean(gt_te_theoretical_per_sample[gt_bool])
    
    ## ======================================================
    ## VT
    ## ======================================================
    tree_type <- ifelse(length(unique(d$Y)) > 2, "regression", "class")
    
    vt.obj <- vt.data(
      dataset         = d_train,
      outcome.field   = "Y",
      treatment.field = "T",
      interactions    = TRUE
    )
    
    vt.for <- vt.forest(
      forest.type  = "one",
      vt.data      = vt.obj,
      interactions = TRUE,
      ntree        = 500
    )
    
    vt.trees <- vt.tree(
      tree.type = tree_type,
      vt.difft  = vt.for,
      threshold = quantile(vt.for$difft, seq(.5, .8, .1)),
      maxdepth  = ncol(d) - 2
    )
    
    vt.sbgrps <- vt.subgroups(vt.trees)
    
    max_te_vt <- -Inf
    max_sg_bool_vt <- rep(TRUE, nrow(d_test))
    best_rule_vt <- "(all)"
    
    if (length(vt.sbgrps$Subgroup) == 0L) {
      # No subgroups: use whole test set
      bool_test <- rep(TRUE, nrow(d_test))
      subgroup_te <- mean(d_test$Y[bool_test & d_test$T == 1]) -
        mean(d_test$Y[bool_test & d_test$T == 0])
      max_te_vt <- subgroup_te
      max_sg_bool_vt <- bool_test
      best_rule_vt <- "(all)"
    } else {
      treatment_effects_vt <- rep(-Inf, length(vt.sbgrps$Subgroup))
      for (i in seq_along(vt.sbgrps$Subgroup)) {
        sg <- vt.sbgrps$Subgroup[i]  # ✅ BUGFIX: use [i], not [1]
        items <- strsplit(sg, " & ")[[1]]
        bool_test <- rep(TRUE, nrow(d_test))
        
        for (item in items) {
          if (grepl(">=", item)) {
            feature_value <- strsplit(item, ">=")[[1]]
            operator <- ">="
          } else {
            feature_value <- strsplit(item, "< ")[[1]]
            operator <- "<"
          }
          feature <- feature_value[1]
          value <- as.numeric(feature_value[2])
          
          if (operator == ">=") {
            bool_test <- bool_test & (d_test[[feature]] >= value)
          } else {
            bool_test <- bool_test & (d_test[[feature]] < value)
          }
        }
        
        subgroup_te <- mean(d_test$Y[bool_test & d_test$T == 1]) -
          mean(d_test$Y[bool_test & d_test$T == 0])
        treatment_effects_vt[i] <- subgroup_te
        
        if (subgroup_te > max_te_vt) {
          max_te_vt <- subgroup_te
          max_sg_bool_vt <- bool_test
          best_rule_vt <- sg
        }
      }
    }
    
    est_te_learned_vt <- mean(d_test$Y[max_sg_bool_vt & d_test$T == 1]) -
      mean(d_test$Y[max_sg_bool_vt & d_test$T == 0])
    gt_te_of_learned_subgroup_vt <-
      mean(gt_te_theoretical_per_sample[max_sg_bool_vt])
    diff_te_GtGt_EstimateLearned_vt <-
      abs(gt_te_of_gt_subgroup - est_te_learned_vt)
    subgroup_size_vt <- sum(max_sg_bool_vt)
    subgroup_hash_vt <- subgroup_hash_fun(max_sg_bool_vt)
    
    ## ======================================================
    ## QUINT
    ## ======================================================
    quint_model <- quint(
      formula = Y ~ T | .,
      data    = d_train
    )
    quint_model_pruned <- prune(quint_model)
    
    if (is.null(quint_model_pruned$si)) {
      subgroups <- NULL
    } else {
      subgroups <- get_rules(quint_model_pruned$si)
    }
    
    max_te_q <- -Inf
    max_sg_bool_q <- rep(TRUE, nrow(d_test))
    best_rule_q <- "(all)"
    
    if (is.null(subgroups)) {
      bool_test <- rep(TRUE, nrow(d_test))
      subgroup_te <- mean(d_test$Y[bool_test & d_test$T == 1]) -
        mean(d_test$Y[bool_test & d_test$T == 0])
      max_te_q <- subgroup_te
      max_sg_bool_q <- bool_test
      best_rule_q <- "(all)"
    } else {
      treatment_effects_q <- rep(-Inf, length(subgroups))
      for (i in seq_along(subgroups)) {
        sg <- subgroups[i]
        items <- strsplit(sg, " & ")[[1]]
        bool_test <- rep(TRUE, nrow(d_test))
        
        for (item in items) {
          if (grepl(">=", item)) {
            feature_value <- strsplit(item, " >= ")[[1]]
            operator <- ">="
          } else {
            feature_value <- strsplit(item, " < ")[[1]]
            operator <- "<"
          }
          feature <- feature_value[1]
          value <- as.numeric(feature_value[2])
          
          if (operator == ">=") {
            bool_test <- bool_test & (d_test[[feature]] >= value)
          } else {
            bool_test <- bool_test & (d_test[[feature]] < value)
          }
        }
        
        subgroup_te <- mean(d_test$Y[bool_test & d_test$T == 1]) -
          mean(d_test$Y[bool_test & d_test$T == 0])
        treatment_effects_q[i] <- subgroup_te
        
        if (subgroup_te > max_te_q) {
          max_te_q <- subgroup_te
          max_sg_bool_q <- bool_test
          best_rule_q <- sg
        }
      }
    }
    
    est_te_learned_q <- mean(d_test$Y[max_sg_bool_q & d_test$T == 1]) -
      mean(d_test$Y[max_sg_bool_q & d_test$T == 0])
    gt_te_of_learned_subgroup_q <-
      mean(gt_te_theoretical_per_sample[max_sg_bool_q])
    diff_te_GtGt_EstimateLearned_q <-
      abs(gt_te_of_gt_subgroup - est_te_learned_q)
    subgroup_size_q <- sum(max_sg_bool_q)
    subgroup_hash_q <- subgroup_hash_fun(max_sg_bool_q)
    
    ## ======================================================
    ## Store results for both methods
    ## ======================================================
    rows_list[[length(rows_list) + 1L]] <- data.frame(
      simulator_name              = simulator_name,
      n                           = n,
      iter                        = iter_,
      method                      = "VT",
      best_rule                   = best_rule_vt,
      est_te_learned              = est_te_learned_vt,
      gt_te                       = gt_te,
      gt_te_of_gt_subgroup        = gt_te_of_gt_subgroup,
      gt_te_of_learned_subgroup   = gt_te_of_learned_subgroup_vt,
      diff_te_GtGt_EstimateLearned = diff_te_GtGt_EstimateLearned_vt,
      subgroup_size               = subgroup_size_vt,
      subgroup_hash               = subgroup_hash_vt,
      stringsAsFactors            = FALSE
    )
    
    rows_list[[length(rows_list) + 1L]] <- data.frame(
      simulator_name              = simulator_name,
      n                           = n,
      iter                        = iter_,
      method                      = "QUINT",
      best_rule                   = best_rule_q,
      est_te_learned              = est_te_learned_q,
      gt_te                       = gt_te,
      gt_te_of_gt_subgroup        = gt_te_of_gt_subgroup,
      gt_te_of_learned_subgroup   = gt_te_of_learned_subgroup_q,
      diff_te_GtGt_EstimateLearned = diff_te_GtGt_EstimateLearned_q,
      subgroup_size               = subgroup_size_q,
      subgroup_hash               = subgroup_hash_q,
      stringsAsFactors            = FALSE
    )
  }
}

sanity_res <- rbindlist(rows_list)

write.csv(sanity_res, "sanity_VT_QUINT_rules.csv", row.names = FALSE)

## ----------------------------------------------------------
## Build a wide comparison VT vs QUINT
## ----------------------------------------------------------

dt <- as.data.table(sanity_res)
dt_wide <- dcast(
  dt,
  simulator_name + n + iter ~ method,
  value.var = c("best_rule",
                "est_te_learned",
                "diff_te_GtGt_EstimateLearned",
                "subgroup_size",
                "subgroup_hash")
)

dt_wide[, same_est_te :=
          est_te_learned_VT == est_te_learned_QUINT]
dt_wide[, same_diff :=
          diff_te_GtGt_EstimateLearned_VT ==
          diff_te_GtGt_EstimateLearned_QUINT]
dt_wide[, same_rule :=
          best_rule_VT == best_rule_QUINT]
dt_wide[, same_hash :=
          subgroup_hash_VT == subgroup_hash_QUINT]

cat("\n=== Sanity check summary (5 iters, n = 3000) ===\n")
print(dt_wide[, .(
  n_pairs           = .N,
  prop_same_est_te  = mean(same_est_te),
  prop_same_diff    = mean(same_diff),
  prop_same_rule    = mean(same_rule),
  prop_same_hash    = mean(same_hash)
)])

cat("\nDetailed per (simulator, iter):\n")
print(dt_wide[, .(
  same_est_te,
  same_diff,
  same_rule,
  same_hash
), by = .(simulator_name, iter)])
