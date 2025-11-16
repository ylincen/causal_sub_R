rm(list=ls())
d = read.csv("./res_cart_new_hidden.csv")
require(data.table)
dt = data.table(d)

# summarize the tree_diff_te
dtt = dt[, .(mean_diff_te = mean(tree_diff_te, na.rm=TRUE),
        sd_diff_te = sd(tree_diff_te, na.rm=TRUE)),
   by=.(n)]
round(dtt$mean_diff_te,3)[seq(1,9,by=2)]
