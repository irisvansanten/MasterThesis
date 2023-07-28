### Description script ###
# Convergence plots were created to check convergence for all scenario's. 
# First a help function is defined, then a function to make the plots and lastly a function to perform the previous function.
# Performed for a model with Z and without Z. Code is the same, only obtaining different file name and saving under different file name. 

### Load libraries ###
library(dplyr)
library(ggplot2)

### Define parameters ###
napproaches = 3
ncat = 4
nconds = 2
nsim = 500
proportions_mv = c(0.1, 0.2, 0.3, 0.4)
true_dist =c (0.1, 0.2, 0.3, 0.4)
MAR_rel = c("strong", "weak")

### Define help functions ###
# 1) Extract average estimated proportions from all simulations
extract_est_proportions = function(prop_boot, nconds, nsim){
  prop_nsim_x = rep(list(vector("list", nsim)), nconds)
  for(cond in 1:nconds){
    for(sim in 1:nsim){
      prop_nsim_x[[cond]][[sim]] = colMeans(bind_rows(prop_boot[[cond]][[sim]]))
    }
  }
  return(prop_nsim_x)
}

########################
### Plot Convergence ###
########################

# Make convergence plots
plot_conv = function(prop_nsim_x_Z, prop_nsim_x_no_Z , MAR_rel, proportions_mv, approach, ncat, nconds, nsim ){
  # Create table with proportion estimate for each simulation
  col_names = c("Condition", "Simulation")
  for(cat in 1:ncat){
    col_names = append(col_names, sprintf("cat_%s", cat))
  }
  df_lv_sim_Z = data.frame(matrix(ncol = 2 + ncat, nrow = nconds*nsim, dimnames = list(NULL, col_names)))
  df_lv_sim_no_Z = data.frame(matrix(ncol = 2 + ncat, nrow = nconds*nsim, dimnames = list(NULL, col_names)))
  count_con = 0
  for(cond in 1:nconds){
    for(sim in 1:nsim){
      count_con = count_con + 1
      df_lv_sim_Z[count_con, 1] = cond
      df_lv_sim_Z[count_con, 2] = sim
      df_lv_sim_no_Z[count_con, 1] = cond
      df_lv_sim_no_Z[count_con, 2] = sim
      count_i = 0
      for(prop in prop_nsim_x_Z[[cond]][[sim]]){
        count_i = count_i + 1
        df_lv_sim_Z[count_con, 2 + count_i] = prop
      }
      count_j = 0
      for(prop2 in prop_nsim_x_no_Z[[cond]][[sim]]){
        count_j = count_j + 1
        df_lv_sim_no_Z[count_con, 2 + count_j] = prop2
      }
    }
    }
  write.csv2(df_lv_sim_Z, file = sprintf("F:/Documents/Thesis_project/lv_sim_approach%s_%s_%s_Z.csv", approach, proportions_mv, MAR_rel), quote = F, row.names = F)
  write.csv2(df_lv_sim_no_Z, file = sprintf("F:/Documents/Thesis_project/lv_sim_approach%s_%s_%s_no_Z.csv", approach, proportions_mv, MAR_rel), quote = F, row.names = F)
  # Create convergence plot 
  for(cond in 1:nconds){
    cond_df_Z = subset(df_lv_sim_Z, Condition == cond)
    cond_df_no_Z = subset(df_lv_sim_no_Z, Condition == cond)
    for(cat in 1:ncat){
      cat_df_Z = cond_df_Z[,c(2, 2 + cat)]
      cat_df_no_Z = cond_df_no_Z[,c(2, 2 + cat)]
      for(row in 2:nrow(cat_df_Z)){
        cat_df_Z[row, 2] = (cat_df_Z[row-1, 2] + cat_df_Z[row, 2])
      }
      for(row2 in 2:nrow(cat_df_no_Z)){
        cat_df_no_Z[row2, 2] = (cat_df_no_Z[row2-1, 2] + cat_df_no_Z[row2, 2])
      }
      cat_df_Z$new_cat = cat_df_Z[, 2]/cat_df_Z[, 1]
      cat_df_no_Z$new_cat = cat_df_no_Z[, 2]/cat_df_no_Z[, 1]
      p = ggplot(data = cat_df_Z, aes(x = Simulation, y = new_cat)) +
        geom_point() +
        theme_minimal()
      ggplot2::ggsave(filename = sprintf("./conv_LC%s_approach%s_%s_%s_cond%s_Z.pdf", cat, approach, proportions_mv, MAR_rel, cond), width = 5, height = 5, plot = p, device = "pdf" )
      p2 = ggplot(data = cat_df_no_Z, aes(x = Simulation, y = new_cat)) +
        geom_point() +
        theme_minimal()
      ggplot2::ggsave(filename = sprintf("./conv_LC%s_approach%s_%s_%s_cond%s_no_Z.pdf", cat, approach, proportions_mv, MAR_rel, cond), width = 5, height = 5, plot = p2, device = "pdf" )
    }
    
  }
  return()
}

# Perform previous defined functions 
check_conv = function(MAR_rel, proportions_mv, napproaches, ncat, nconds, nsim){
  for(i in MAR_rel){
    for(prop_mv in proportions_mv){
      for(approach in 1:napproaches){
        prop_boot_Z = readRDS(sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_Z.rds", prop_mv, i))
        prop_nsim_x_Z = extract_est_proportions(prop_boot_Z[[approach]], nconds, nsim)
        prop_boot_no_Z = readRDS(sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_no_Z.rds", prop_mv, i))
        prop_nsim_x_no_Z = extract_est_proportions(prop_boot_no_Z[[approach]], nconds, nsim)
        plot_conv(prop_nsim_x_Z, prop_nsim_x_no_Z, i, prop_mv, approach, ncat, nconds, nsim)
      }
    }
  }
}
check_conv(MAR_rel, proportions_mv, napproaches, ncat, nconds, nsim)
