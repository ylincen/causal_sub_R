rm(list=ls())

require(data.table)
require(dplyr)

d_cart = fread("./semi_cart_results.csv")
d_rsides = fread("./semi_rsides_results.csv")
options(warn=1)

library(ggplot2)

#—Option 1: direct (closest to your original) ----
ggplot() +
  geom_point(
    aes(x = d_cart$max_subgroup_treatment_effect,
        y = d_cart$max_subgroup_gt_treatment_effct,
        colour = "CART"),
    size = 2
  ) +
  geom_point(
    aes(x = d_rsides$max_subgroup_treatment_effect,
        y = d_rsides$max_subgroup_gt_treatment_effct,
        colour = "R-SIDES"),
    size = 2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color="red") +
  scale_colour_manual(values = c("CART" = "black", "R-SIDES" = "blue"), name = NULL) +
  labs(x = "Estimated subgroup treatment effect", y = "Ground-truth subgroup treatment effect") +
  theme_bw()

if(nrow(d_cart) != nrow(d_rsides)){
  d_cart = d_cart[1:nrow(d_rsides),]
}
plot(d_cart$max_subgroup_treatment_effect, d_cart$max_subgroup_gt_treatment_effct)
abline(a=0, b=1, col="red")
points(d_rsides$max_subgroup_treatment_effect, d_rsides$max_subgroup_gt_treatment_effct, col="blue")

mean(d_cart$max_subgroup_treatment_effect > d_rsides$max_subgroup_treatment_effect)

cat("mean/sd of d_cart max_subgroup_treatment_effect: ", 
      mean(d_cart$max_subgroup_treatment_effect, na.rm = TRUE), 
      sd(d_cart$max_subgroup_treatment_effect, na.rm = TRUE))

cat("mean/sd of d_rsides max_subgroup_treatment_effect: ",
      mean(d_rsides$max_subgroup_treatment_effect, na.rm = TRUE), 
      sd(d_rsides$max_subgroup_treatment_effect, na.rm = TRUE))

cat("mean/sd of d_cart max_subgroup_diff_average: ",
      mean(d_cart$max_subgroup_diff_average, na.rm = TRUE), 
      sd(d_cart$max_subgroup_diff_average, na.rm = TRUE))
cat("mean/sd of d_rsides max_subgroup_diff_average: ",
      mean(d_rsides$max_subgroup_diff_average, na.rm = TRUE), 
      sd(d_rsides$max_subgroup_diff_average, na.rm = TRUE))

compare = cbind(d_cart[,c("dataset", "max_subgroup_treatment_effect", 
                "max_subgroup_gt_treatment_effct")],
      d_rsides[,c("dataset", "max_subgroup_treatment_effect", 
                  "max_subgroup_gt_treatment_effct")])
colnames(compare) = c("dataset", 
                            "cart_max_subgroup_treatment_effect", 
                            "cart_max_subgroup_gt_treatment_effct",
                      "dataset2",
                            "rsides_max_subgroup_treatment_effect", 
                            "rsides_max_subgroup_gt_treatment_effct")
all(compare$dataset == compare$dataset2)
compare$dataset2 = NULL
compare$cart_win = compare$cart_max_subgroup_treatment_effect > 
  compare$rsides_max_subgroup_treatment_effect


# plot(d_cart$max_subgroup_treatment_effect, 
#      d_rsides$max_subgroup_treatment_effect)
# plot(d_cart$max_subgroup_diff_average, 
#      d_rsides$max_subgroup_diff_average)
