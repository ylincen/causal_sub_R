rm(list=ls())
require(aVirtualTwins)
d = read.csv("./new_simulation/simulate1/n_1000_iter_2.csv")

train_indices = sample(1:nrow(d), size = 0.5 * nrow(d), replace = FALSE)
d_train = d[train_indices, ]
d_test = d[-train_indices, ]


vt.obj <- vt.data(dataset         = d_train,
                  outcome.field   = "Y",
                  treatment.field = "T",
                  interactions    = TRUE)
# Print Incidences of sepsis data
vt.obj$getIncidences()

# First step : create random forest model
vt.for <- vt.forest(forest.type  = "one",
                    vt.data      = vt.obj,
                    interactions = TRUE,
                    ntree        = 500)
# Second step : find rules in data 
vt.trees <- vt.tree(tree.type = "class",
                    vt.difft  = vt.for, 
                    threshold = quantile(vt.for$difft, seq(.5,.8,.1)),
                    maxdepth  = 5)
# Print results
vt.sbgrps <- vt.subgroups(vt.trees)

# get the subgroup analysis results
max_te = -Inf
gt_bool = d_test$X1 > 1
gt_te = mean(d_test$Y[gt_bool & (d_test$T == 1)]) - 
  mean(d_test$Y[gt_bool & (d_test$T == 0)])


for(i in 1:length(vt.sbgrps$Subgroup)){
  sg = vt.sbgrps$Subgroup[1]
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
  if(subgroup_te > max_te){
    max_te = subgroup_te
    max_subgroup_jaccard = sum(bool_test & gt_bool) / 
      sum(bool_test | gt_bool)
    abs_error_te = abs(max_te - gt_te)
  }
}

cat("Maximum treatment effect in subgroup analysis: ", max_te, "\n",
    "Jaccard index with ground truth subgroup: ", max_subgroup_jaccard, "\n",
    "Absolute error of treatment effect: ", abs_error_te, "\n")

