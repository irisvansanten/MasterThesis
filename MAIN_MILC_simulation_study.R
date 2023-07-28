### Description script ###
# Main script to perform simulation study. 
# First some preparation: empty environment, load libraries, define parameters and define population values
# The main code is divided in 1) data generation, 2) introduce and treat missing values, 
# 3) perform MILC, 4) obtain results.
# Another script was loaded in parallelizing part of the code, FUN_help.R 

### Empty environment ###
rm(list = ls())
setwd('F://Documents//Thesis_project')

### Load libraries ###
library(Rlab)
library(poLCA)
library(confreq)
library(dplyr)
library(tidyr)
library(broom)
library(resample)
library(ggplot2)
library(gridExtra)
library(missForest)
library(mice)
library(tidyverse)
library(scales)
library(parallel)
library(doParallel)
library(doRNG)
library(bigstatsr)
source("FUN_help.R")

### Define parameters ###
napproaches = 3
nboot = 5
ncat = 4
nconds = 2
nsim = 500
nsize = 10000
MAR_rel_prop = list(c(0.2, 0.5, 0.01, 0.6), c(0.4, 0.3, 0.35, 0.3))
proportions_mv = c(0.1, 0.2, 0.3, 0.4)
T_dist = c(0.1, 0.2, 0.3, 0.4)
lv_pm = c("lv", "RMSE", "var", "bias")
MAR_rel = c("strong", "weak")
ME_levels =  c("5%", "20%")

### Define T population probabilities ###
# Probability distributions of T population that has either 5% measurement error or 20% measurement error. Data is simulated based on these matrices. Each matrix refers to a variable. 
# First four identical matrices refer to Y1 to Y4 and the last refers to Z. 
# Each matrix contains 4 rows and columns representing the probability to observe a certain category (rows) given the T latent classes (columns). 
pop_prob = list(
  # 5% ME
  list(matrix(c(285/300, 5/300, 5/300, 5/300,
                5/300, 285/300, 5/300, 5/300, 
                5/300, 5/300, 285/300, 5/300,
                5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow =  T),
       matrix(c(285/300, 5/300, 5/300, 5/300,
                5/300, 285/300, 5/300, 5/300, 
                5/300, 5/300, 285/300, 5/300,
                5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow =  T),
       matrix(c(285/300, 5/300, 5/300, 5/300,
                5/300, 285/300, 5/300, 5/300, 
                5/300, 5/300, 285/300, 5/300,
                5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow =  T),
       matrix(c(285/300, 5/300, 5/300, 5/300,
                5/300, 285/300, 5/300, 5/300, 
                5/300, 5/300, 285/300, 5/300,
                5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow =  T),
       matrix(c(18/30, 4/30, 4/30, 4/30,
                8/30, 6/30, 8/30, 8/30, 
                2/30, 2/30, 24/30, 2/30,
                5/30, 5/30, 5/30, 15/30), ncol = ncat, byrow =  T)
  ), 
  #20%
  list(matrix(c(24/30, 2/30, 2/30, 2/30,
                2/30, 24/30, 2/30, 2/30,
                2/30, 2/30, 24/30, 2/30,
                2/30, 2/30, 2/30, 24/30), ncol = ncat, byrow =  T),
       matrix(c(24/30, 2/30, 2/30, 2/30,
                2/30, 24/30, 2/30, 2/30,
                2/30, 2/30, 24/30, 2/30,
                2/30, 2/30, 2/30, 24/30), ncol = ncat, byrow =  T),
       matrix(c(24/30, 2/30, 2/30, 2/30,
                2/30, 24/30, 2/30, 2/30,
                2/30, 2/30, 24/30, 2/30,
                2/30, 2/30, 2/30, 24/30), ncol = ncat, byrow =  T),
       matrix(c(24/30, 2/30, 2/30, 2/30,
                2/30, 24/30, 2/30, 2/30,
                2/30, 2/30, 24/30, 2/30,
                2/30, 2/30, 2/30, 24/30), ncol = ncat, byrow =  T),
       matrix(c(18/30, 4/30, 4/30, 4/30,
                8/30, 6/30, 8/30, 8/30, 
                2/30, 2/30, 24/30, 2/30,
                5/30, 5/30, 5/30, 15/30), ncol = ncat, byrow =  T)
  )
)

#######################
#######################
### Data generation ###
#######################
#######################

### Main function ###
data_generation = function(T_dist, ncat, nconds, nsim, nsize, pop_prob){
  simboot_dat = vector("list", nconds)
  model_mv = vector("list", nconds)
  for(cond in 1:nconds){
    for(sim in 1:nsim){
      set.seed(17082023)
      model = poLCA.simdata(N       = nsize,
                            nclass  = ncat,
                            probs   = pop_prob[[cond]],
                            P       = T_dist,
                            missval = F)
      model_mv[[cond]][[sim]] = model
      # Change from numeric to categoric
      model_mv[[cond]][[sim]]$dat$Y1 = factor(model_mv[[cond]][[sim]]$dat$Y1, levels = 1:5)
      model_mv[[cond]][[sim]]$dat$Y2 = factor(model_mv[[cond]][[sim]]$dat$Y2, levels = 1:5)
      model_mv[[cond]][[sim]]$dat$Y3 = factor(model_mv[[cond]][[sim]]$dat$Y3, levels = 1:5)
      model_mv[[cond]][[sim]]$dat$Y4 = factor(model_mv[[cond]][[sim]]$dat$Y4, levels = 1:5)
      model_mv[[cond]][[sim]]$dat$Y5 = factor(model_mv[[cond]][[sim]]$dat$Y5, levels = 1:5)
    }
  }
  saveRDS(model_mv, file = sprintf("F:/Documents/Thesis_project/model_mv.rds"))
  return(model_mv)
}

### Perform functions ###
model_mv = data_generation(T_dist, ncat, nconds, nsim, nsize, pop_prob)

##########################################
##########################################
### Introduce and treat missing values ###
##########################################
##########################################

### Main function ###
missing_data = function(model, MAR_rel, MAR_rel_prop, proportions_mv, napproaches, nboot, nconds, nsim){
  if(MAR_rel == "weak"){
    cat_miss_prop = MAR_rel_prop[[2]]
  }else if(MAR_rel == "strong"){
    cat_miss_prop = MAR_rel_prop[[1]]
  }
  starttime = Sys.time()
  nCores = detectCores()-1
  myCluster = parallel::makeCluster(nCores, type = "PSOCK")
  registerDoParallel(myCluster)
  clusterSetRNGStream(cl = myCluster, iseed = 666)
  clusterExport(cl = myCluster, varlist=c('ncat', 'napproaches', 'nboot', 'nconds', 'nsim'))
  simboot_dat_mv_a = vector("list", nconds)
  simboot_dat_mv_b = vector("list", nconds)
  for(cond in 1:nconds){
    simboot_dat_mv_parallel = foreach(sim = 1:nsim,
                                      .packages = c("poLCA", "confreq", "dplyr", "resample", "parallel", "bigstatsr", "Rlab", 'mice') )%dorng%{
                                        source("FUN_help.R")
                                        simboot_dat_a =  vector("list", napproaches)
                                        simboot_dat_b =  vector("list", napproaches)
                                        dataframe = missing_values(model[[cond]][[sim]],
                                                                   cat_miss_prop,
                                                                   proportions_mv,
                                                                   ncat)
                                        for(approach in 1:napproaches){
                                          dataframe_mv = approaches_function(dataframe$dataframe_NA,
                                                                             dataframe$dataframe_model_NA,
                                                                             dataframe$model,
                                                                             approach)
                                          dffreq_mv_a = dataframe_mv[, 1:ncol(dataframe_mv)] %>%
                                            dplyr::count(Y1, Y2, Y3, Y4, Y5,  Tclass)
                                          simboot_dat_a[[approach]] = bootdata(data = dffreq_mv_a, nboot)
                                          dffreq_mv_b = dataframe_mv[, 1:ncol(dataframe_mv)] %>%
                                            dplyr::count(Y1, Y2, Y3, Y4,  Tclass)
                                          simboot_dat_b[[approach]] = bootdata(data = dffreq_mv_b, nboot)
                                        }
                                        return(list(a = simboot_dat_a, b = simboot_dat_b) )
                                      }
    a_part = rep(list(vector("list", napproaches)), nsim)
    b_part = rep(list(vector("list", napproaches)), nsim)
    for(sim in 1:nsim){
      for(approach in 1:napproaches){
        a_part[[sim]][[approach]] = simboot_dat_mv_parallel[[sim]]$a[[approach]]
        b_part[[sim]][[approach]] = simboot_dat_mv_parallel[[sim]]$b[[approach]]
      }
      
    }
    simboot_dat_mv_a[[cond]] = a_part
    simboot_dat_mv_b[[cond]] = b_part
  }
  stopCluster(myCluster)
  SimulationTotalTime = Sys.time() - starttime
  print(SimulationTotalTime)
  saveRDS(simboot_dat_mv_a, file = sprintf("F:/Documents/Thesis_project/simboot_dat_mv_%s_%s_Z.rds", proportions_mv, MAR_rel))
  saveRDS(simboot_dat_mv_b, file = sprintf("F:/Documents/Thesis_project/simboot_dat_mv_%s_%s_no_Z.rds", proportions_mv, MAR_rel))
  return()
}

### Perform functions ###
model_mv = readRDS("~/Thesis_project/model_mv.rds")
MV_approach = function(model_mv, MAR_rel, MAR_rel_prop, proportions_mv, napproaches, nboot, nconds, nsim){
  lapply(MAR_rel, function(MAR_rel){
    lapply(proportions_mv,
           function(proportions_mv){missing_data(model_mv, MAR_rel, MAR_rel_prop, proportions_mv, napproaches, nboot, nconds, nsim)}
    )
  })
}
dataset_MV = MV_approach(model_mv, MAR_rel, MAR_rel_prop, proportions_mv, napproaches, nboot, nconds, nsim)

##################################
##################################
### Perform MILC for 2 models ####
##################################
##################################

### Main function for model with Z ###
MILC_fuction_no_Z = function(simbootdat_mv, MAR_rel, proportions_mv, napproaches, nboot, nconds, nsim){
  starttime = Sys.time()
  nCores = detectCores() - 1
  myCluster = parallel::makeCluster(nCores, type = "PSOCK")
  registerDoParallel(myCluster)
  clusterSetRNGStream(cl = myCluster, iseed = 666)
  clusterExport(cl = myCluster, varlist = c('ncat', 'napproaches', 'nboot', 'nconds', 'nsim'))
  prop_boot = vector("list", napproaches)
  for(approach in 1:napproaches){
    prop_boot_2 = vector("list", nconds)
    for(cond in 1:nconds){
      result_MILC = foreach(sim = 1:nsim, 
                            .packages=c("poLCA", "confreq", "dplyr", "resample", "parallel", "bigstatsr", "Rlab", 'mice')) %dorng%{
                              source("FUN_help.R")
                              impdats = vector("list", nboot) 
                              prop_boot_mv = vector("list", nboot) 
                              for(boot in 1:nboot){
                                log = capture.output({
                                  impdats[[boot]] = posteriors_no_Z(bootdat = simbootdat_mv[[cond]][[sim]][[approach]], k = boot)})
                                for (l in 1:(sum(simbootdat_mv[[cond]][[sim]][[approach]]$n))){
                                  impdats[[boot]][l, "imp"] = which(rmultinom(n = 1, size = 1,
                                                                              prob = impdats[[boot]][l,c("1","2","3","4")]) == 1)
                                }
                                prop_boot_mv[[boot]] =
                                  prop.table(table(impdats[[boot]]$imp))
                              }
                              return(prop_boot_mv)
                            }
      prop_boot_2[[cond]] = result_MILC
    }
    prop_boot[[approach]] = prop_boot_2
  }
  stopCluster(myCluster)
  SimulationTotalTime = Sys.time() - starttime
  print(SimulationTotalTime)
  saveRDS(prop_boot, file = sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_no_Z.rds", proportions_mv, MAR_rel))
  return(prop_boot)
}

### Main function for model without Z ###
MILC_fuction_Z = function(simbootdat_mv, MAR_rel, proportions_mv, napproaches, nboot, nconds, nsim){
  starttime = Sys.time()
  nCores = detectCores() - 1
  myCluster = parallel::makeCluster(nCores, type="PSOCK")
  registerDoParallel(myCluster)
  clusterSetRNGStream(cl = myCluster, iseed = 666)
  clusterExport(cl = myCluster, varlist=c('ncat', 'napproaches', 'nboot', 'nconds', 'nsim'))
  prop_boot = vector("list", napproaches)
  for(approach in 1:napproaches){
    prop_boot_2 = vector("list", nconds)
    for(cond in 1:nconds){
      result_MILC = foreach(sim = 1:nsim, 
                            .packages = c("poLCA", "confreq", "dplyr", "resample", "parallel", "bigstatsr", "Rlab", 'mice')) %dorng%{
                              source("FUN_help.R")
                              impdats = vector("list", nboot) 
                              prop_boot_mv = vector("list", nboot) 
                              for(boot in 1:nboot){
                                log = capture.output({
                                  impdats[[boot]] = posteriors_Z(bootdat = simbootdat_mv[[cond]][[sim]][[approach]], k = boot)})
                                for (l in 1:(sum(simbootdat_mv[[cond]][[sim]][[approach]]$n))){
                                  impdats[[boot]][l, "imp"] = which(rmultinom(n = 1, size = 1,
                                                                              prob = impdats[[boot]][l,c("1","2","3","4")]) == 1)
                                }
                                prop_boot_mv[[boot]] =
                                  prop.table(table(impdats[[boot]]$imp))
                              }
                              return(prop_boot_mv)
                            }
      prop_boot_2[[cond]] = result_MILC
    }
    prop_boot[[approach]] = prop_boot_2
  }
  stopCluster(myCluster)
  SimulationTotalTime = Sys.time() - starttime
  print(SimulationTotalTime)
  saveRDS(prop_boot, file = sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_Z.rds", proportions_mv, MAR_rel))
  return(prop_boot)
}

### Perform functions ###
MILC_MV = function(MAR_rel, proportions_mv, napproaches, nboot, nconds, nsim){
  for(i in MAR_rel){
    for(prop_mv in proportions_mv){
      simbootdat_mv_no_Z = readRDS(sprintf("~/Thesis_project/simboot_dat_mv_%s_%s_no_Z.rds", prop_mv, i))
      prop_boot_no_Z = MILC_fuction_no_Z(simbootdat_mv_no_Z, i, prop_mv, napproaches, nboot, nconds, nsim)
      simbootdat_mv_Z = readRDS(sprintf("~/Thesis_project/simboot_dat_mv_%s_%s_Z.rds", prop_mv, i))
      prop_boot_Z = MILC_fuction_Z(simbootdat_mv_Z, i, prop_mv, napproaches, nboot, nconds, nsim)
    }
  }
}
MILC = MILC_MV(MAR_rel, proportions_mv, napproaches, nboot, nconds, nsim)

#######################
#######################
### Obtain results ####
#######################
#######################

### Help function ###
# 1) Obtain pooled estimated latent class proportions of each bootstrap sample
extract_est_proportions = function(prop_boot, nconds, nsim){
  prop_nsim_x = rep(list(vector("list", nsim)), nconds)
  for(cond in 1:nconds){
    for(sim in 1:nsim){
      prop_nsim_x[[cond]][[sim]] = colMeans(bind_rows(prop_boot[[cond]][[sim]]))
    }
  }
  return(prop_nsim_x)
}

# 2) obtain average estimates over all simulations
lv_pm_result = function(df, col_names, ncat, nconds){
  new_df = data.frame(matrix(ncol = 1 + ncat, nrow = nconds, dimnames = list(NULL, col_names)))
  count_con = 0
  for(cond in 1:nconds){
    count_con = count_con + 1
    new_df[count_con, 1] = cond
    for(col in 2:ncol(new_df)){
      new_df[count_con, col] = mean(df[df$Condition == cond, ][, col+1])
    }
  }
  return(new_df)
}

# 3) Create table with estimates 
est_of_interest = function(prop_nsim_x, ncat, nconds, nsim ){
  # Table with estimate of interest for each simulation
  col_names_sim_lv = c("Condition", "Simulation")
  for(cat in 1:ncat){
    col_names_sim_lv = append(col_names_sim_lv, sprintf("cat_%s", cat))
  }
  df_lv_sim = data.frame(matrix(ncol = 2 + ncat, nrow = nconds*nsim, dimnames = list(NULL, col_names_sim_lv)))
  count_con = 0
  for(cond in 1:nconds){
    for(sim in 1:nsim){
      count_con=count_con + 1
      df_lv_sim[count_con, 1] = cond
      df_lv_sim[count_con,2] = sim
      count_i = 0
      for(prop in prop_nsim_x[[cond]][[sim]]){
        count_i = count_i+1
        df_lv_sim[count_con, 2 + count_i] = prop
      }
    }
  }
  col_names_lv = c("Condition")
  for(cat in 1:ncat){
    col_names_lv = append(col_names_lv, sprintf("cat_%s",cat))
  }
  table_lv = lv_pm_result(df_lv_sim, col_names_lv, ncat, nconds)
  return(list("df_lv_sim" = df_lv_sim, "table_lv" = table_lv))
}

# 4) Calculate the performance measures 
performance_measures = function(df_lv_sim, table_lv, MAR_rel, proportions_mv, T_dist, approach, ncat, nconds, nsim){
  # Table bias: per simulation and pooled 
  col_names_sim = c("Condition", "Simulation")
  for(cat in 1:ncat){
    col_names_sim = append(col_names_sim, sprintf("cat_%s", cat))
  }
  df_bias_sim = data.frame(matrix(ncol = 2+ncat, nrow = nsim*nconds, dimnames = list(NULL, c(col_names_sim))))
  df_bias_sim$Condition = df_lv_sim$Condition
  df_bias_sim$Simulation = df_lv_sim$Simulation
  for(cat in 1:ncat){
    df_bias_sim[, 2 + cat] = df_lv_sim[, 2 + cat] - T_dist[cat] 
  }
  col_names = c("Condition")
  for(cat in 1:ncat){
    col_names = append(col_names, sprintf("cat_%s", cat))
  }
  table_bias = lv_pm_result(df_bias_sim, col_names, ncat, nconds)
  # Table RMSE
  ## MSE
  df_MSE_sim = data.frame(matrix(ncol = 2+ncat, nrow = nsim*nconds, dimnames = list(NULL, c(col_names_sim))))
  df_MSE_sim$Condition = df_lv_sim$Condition
  df_MSE_sim$Simulation = df_lv_sim$Simulation
  df_MSE_sim[, 3:(ncat+2)] = (df_bias_sim[, 3:(ncat+2)])^2 
  table_MSE = lv_pm_result(df_MSE_sim, col_names, ncat, nconds)
  ## RMSE
  table_RMSE = data.frame(matrix(ncol = 1+ncat, nrow = nconds, dimnames = list(NULL, c(col_names))))
  table_RMSE$Condition = table_lv$Condition
  RMSE = apply(table_MSE[, 2:(ncat+1)], 2, sqrt) #square root of MSE
  table_RMSE[, 2:(ncat+1)] = RMSE
  
  # Table variance
  table_var = data.frame(matrix(ncol = 1+ncat, nrow = nconds, dimnames = list(NULL, c(col_names))))
  table_var$Condition = table_lv$Condition
  var = (table_RMSE[, 2:(ncat+1)])^2 - (table_bias[, 2:(ncat+1)])^2 # var = MSE -  Bias^2 
  #var = table_MSE[,2:(ncat+1)] - (table_bias[,2:(ncat+1)])^2 --> gives same result
  table_var[, 2:(ncat+1)] = var
  return(list("table_bias" = table_bias, "table_RMSE" = table_RMSE, "table_var" = table_var))
}

### Perform functions ###
MILC_results = function(MAR_rel, proportions_mv, T_dist, napproaches, ncat, nconds, nsim){
  for(i in MAR_rel){
    for(prop_mv in proportions_mv){
      for(approach in 1:napproaches){
        prop_boot_no_Z = readRDS(sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_no_Z.rds", prop_mv, i))
        prop_nsim_x_no_Z = extract_est_proportions(prop_boot_no_Z[[approach]], nconds, nsim)
        lv_results_no_Z = est_of_interest(prop_nsim_x_no_Z, ncat, nconds, nsim)
        write.csv2(lv_results_no_Z$table_lv, file = sprintf("F:/Documents/Thesis_project/table_lv_approach%s_%s_%s_no_Z.csv", approach, prop_mv, i), quote=F, row.names=F)
        pm_results_no_Z = performance_measures(lv_results_no_Z$df_lv_sim, lv_results_no_Z$table_lv, i, prop_mv, T_dist, approach, ncat, nconds, nsim)
        write.csv2(pm_results_no_Z$table_bias, file = sprintf("F:/Documents/Thesis_project/table_bias_approach%s_%s_%s_no_Z.csv", approach, prop_mv, i), quote=F, row.names=F ) 
        write.csv2(pm_results_no_Z$table_RMSE, File = sprintF("F:/Documents/Thesis_project/table_RMSE_approach%s_%s_%s_no_Z.csv", approach, prop_mv, i), quote=F,row.names=F )
        write.csv2(pm_results_no_Z$table_var, file = sprintf("F:/Documents/Thesis_project/table_var_approach%s_%s_%s_no_Z.csv", approach, prop_mv, i), quote=F, row.names=F )
        
        prop_boot_Z = readRDS(sprintf("F:/Documents/Thesis_project/prop_boot_approach_%s_%s_Z.rds", prop_mv, i))
        prop_nsim_x_Z = extract_est_proportions(prop_boot_Z[[approach]], nconds, nsim)
        lv_results_Z = est_of_interest(prop_nsim_x_Z, ncat, nconds, nsim)
        write.csv2(lv_results_Z$table_lv, file=sprintf("F:/Documents/Thesis_project/table_lv_approach%s_%s_%s_Z.csv", approach, prop_mv, i), quote=F, row.names=F)
        pm_results_Z = performance_measures(lv_results_Z$df_lv_sim, lv_results_Z$table_lv, i, prop_mv, T_dist, approach, ncat, nconds, nsim)
        write.csv2(pm_results_Z$table_bias, file = sprintf("F:/Documents/Thesis_project/table_bias_approach%s_%s_%s_Z.csv", approach, prop_mv, i), quote=F, row.names=F ) 
        write.csv2(pm_results_Z$table_RMSE, file = sprintf("F:/Documents/Thesis_project/table_RMSE_approach%s_%s_%s_Z.csv", approach, prop_mv, i), quote=F,row.names=F )
        write.csv2(pm_results_Z$table_var, file = sprintf("F:/Documents/Thesis_project/table_var_approach%s_%s_%s_Z.csv", approach, prop_mv, i), quote=F, row.names=F )
      }
    }
  }
}
dataset_MILC = MILC_results(MAR_rel, proportions_mv, T_dist, napproaches, ncat, nconds, nsim)

#####################
#####################
### Plot results ####
#####################
#####################

### Help function ###
# 1) Function to calculate max and minimum value to set the y limits the same in all conditions
min_and_max = function(name){
  mylist = list.files("F:/Documents/Thesis_project/", pattern = sprintf("%s", name))
  df_combined = list()
  for(i in mylist){
    table = read.csv(file = sprintf("Tables\\%s", i ), header = T, sep = ";")
    table_est = table[, grep("cat", colnames(table))]
    df_combined[[i]] = table_est
  }
  df_final = do.call(rbind, df_combined )
  df_final[1:4] = lapply(df_final[1:4], gsub, pattern = ",", replacement = ".")
  df_final[] = as.numeric(as.matrix(df_final))
  if(name == "lv"){
    my_list = list("min1" = min(df_final[, 1]), "min2" = min(df_final[, 2]), "min3" = min(df_final[, 3]), "min4" = min(df_final[, 4]),
                   "max1" = max(df_final[, 1]), "max2" = max(df_final[, 2]), "max3" = max(df_final[, 3]), "max4" = max(df_final[, 4]))
  }else{
    my_list = list("min" = min(df_final), "max" = max(df_final))
  }
  return(my_list)
}

### Main function ###
plotting_results = function(ME_levels, cat, rel, proportions_mv, napproaches, nconds, pm_name = 0){
  min_max = min_and_max(pm_name)
  if(pm_name == "lv"){
    max_var = min_and_max("var")[[2]]
    max_sd = sqrt(max_var)*1.96
    y_limit = c(min_max[[cat]] - max_sd,min_max[[cat+4]] + max_sd)
  }else{
    y_limit = c(min_max$min, min_max$max)
  }
  for(approach in 1:napproaches){
    models = c("no_Z", "Z")
    ME_condition = c()
    y_axis = c()
    var = c()
    model = c(rep("no_Z", times = length(proportions_mv)*length(ME_levels)),rep("Z", times = length(proportions_mv)*length(ME_levels)) )
    for(me in 1: length(ME_levels)){
      for(x in proportions_mv){
        table = read.csv(file = sprintf("Tables\\table_%s_approach%s_%s_%s_no_Z.csv", pm_name, approach, x, rel ), header=T, sep=";")
        table_est = table[, grep("cat", colnames(table))]
        y_axis = append(y_axis, as.numeric(sub(",", ".", table_est[me, cat], fixed = T)))
        
        table_Z = read.csv(file = sprintf("Tables\\table_%s_approach%s_%s_%s_Z.csv", pm_name, approach, x, rel ), header=T, sep=";")
        table_est_Z = table_Z[, grep("cat", colnames(table_Z))]
        y_axis = append(y_axis, as.numeric(sub(",", ".", table_est_Z[me, cat], fixed = T)))
        
        table_var = read.csv(file = sprintf("Tables\\table_var_approach%s_%s_%s_no_Z.csv", approach, x, rel ), header=T, sep=";")
        table_var_est = table_var[, grep("cat", colnames(table_var))]
        var = append(var, as.numeric(sub(",", ".", table_var_est[me, cat], fixed = T)))
        
        table_var_Z = read.csv(file = sprintf("Tables\\table_var_approach%s_%s_%s_Z.csv", approach, x, rel ), header=T, sep=";")
        table_var_est_Z = table_var_Z[, grep("cat", colnames(table_var_Z))]
        var = append(var, as.numeric(sub(",", ".", table_var_est_Z[me, cat], fixed = T)))
        ME_condition = append(ME_condition, ME_levels[me])
      }}
    comb = paste(model, factor(ME_condition, levels = ME_levels), sep = "")
    plot = data.frame("x" = rep(proportions_mv, nconds), "model" = model, 
                      "ME_levels" = factor(ME_condition, levels = ME_levels), 
                      "y" = y_axis, "sd" = sqrt(var)*1.96, "comb" = comb)
    write.csv2(plot, file = sprintf("F:/Documents/Thesis_project/table_%s_LC%s_%s_approach%s_plots.csv", pm_name,  cat, rel, approach), quote = F, row.names = F )
    if(pm_name == "lv"){
      # Results for appendix
      p_appendix = ggplot(plot, aes(x = x*100, y = y, colour = comb, shape = comb, group = interaction(model, ME_levels))) + 
        geom_point() + 
        geom_errorbar(aes(ymin = y - sd, ymax = y + sd), alpha = 0.5) + 
        geom_line(linetype = "dashed", alpha = 0.4) +
        ylim(y_limit) + 
        xlab("Missing data (%)") + 
        ylab("Estimates") +
        theme_minimal() +
        scale_colour_manual(name = "Model and ME",
                            labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                            values = c("blue4", "darkred", "blue", "red")) +   
        scale_shape_manual(name = "Model and ME",
                           labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                           values = c(15, 15, 17,17))+
        theme(text = element_text(size = 14), legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2, byrow = T))+
        ggtitle(sprintf("")) 
      ggplot2::ggsave(filename = sprintf("./%s_LC%s_%s_approach%s_appendix.pdf", pm_name,  cat, rel, approach), width = 5, height = 5, plot = p_appendix, device = "pdf" )
      # Results for in main part report
      p_main = ggplot(plot, aes(x = x*100, y = y, colour = comb, shape = comb,
                                group = interaction(model, ME_levels))) + 
        geom_point() + 
        geom_errorbar(aes(ymin = y - sd, ymax = y + sd), alpha = 0.5) + 
        geom_line(linetype = "dashed", alpha = 0.4)+
        xlab("Missing data (%)") + 
        ylab("Estimates" ) +
        geom_hline(yintercept = cat/10)+
        theme_minimal() +
        scale_colour_manual(name = "Model and ME",
                            labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                            values = c("blue4", "darkred", "blue", "red")) +   
        scale_shape_manual(name = "Model and ME",
                           labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                           values = c(15, 15, 17,17))+
        theme(text = element_text(size = 14), legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2, byrow = T))+
        ggtitle(sprintf("")) 
      ggplot2::ggsave(filename = sprintf("./%s_LC%s_%s_approach%s_main.pdf", pm_name,  cat, rel, approach), width = 5, height = 5, plot = p_main, device = "pdf" )
    }else{
      if(pm_name == "bias"){
        pm_name = "Bias"
      }
      # Results for appendix
      p_appendix = ggplot(plot, aes(x = x*100, y = y, colour = comb, shape = comb,
                                    group = interaction(model, ME_levels))) + 
        geom_point() + 
        geom_line(linetype = "dashed", alpha = 0.4)+
        ylim(y_limit) + 
        xlab("Missing data (%)") + 
        ylab(sprintf("%s", pm_name)) +
        theme_minimal() +
        scale_colour_manual(name = "Model and ME",
                            labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                            values = c("blue4", "darkred", "blue", "red")) +   
        scale_shape_manual(name = "Model and ME",
                           labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                           values = c(15, 15, 17, 17))+
        theme(text = element_text(size = 14), legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2, byrow = T)) +
        ggtitle(sprintf("")) 
      if(pm_name == "Bias"){
        pm_name = "bias"
      }
      ggplot2::ggsave(filename = sprintf("./%s_LC%s_%s_approach%s_appendix.pdf", pm_name,  cat, rel, approach), width = 5, height = 5, plot = p_appendix, device = "pdf" )
      # Results for in main part report
      if(pm_name == "bias"){
        pm_name = "Bias"
      }
      p_main = ggplot(plot, aes(x = x*100, y = y, colour = comb, shape = comb,
                                group = interaction(model, ME_levels))) + 
        geom_point() + 
        geom_line(linetype = "dashed", alpha = 0.4)+
        xlab("Missing data (%)") + 
        ylab(sprintf("%s", pm_name) ) +
        theme_minimal() +
        scale_colour_manual(name = "Model and ME",
                            labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                            values = c("blue4", "darkred", "blue", "red")) +   
        scale_shape_manual(name = "Model and ME",
                           labels = c("no Z, 20%", "no Z, 5%", "Z, 20%", "Z, 5%"),
                           values = c(15, 15, 17,17))+
        theme(text = element_text(size = 14), legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2, byrow = T))+
        ggtitle(sprintf(""))
      if(pm_name == "Bias"){
        pm_name = "bias"
      }
      ggplot2::ggsave(filename = sprintf("./%s_LC%s_%s_approach%s_main.pdf", pm_name,  cat, rel, approach), width = 5, height = 5, plot = p_main, device = "pdf" )
    }
  }
  
}

### Perform functions ###
results = function(ME_levels, lv_pm, MAR_rel, proportions_mv, napproaches, ncat, nconds){
  for(rel in 1:length(MAR_rel)){
    for(cat in 1:ncat){
      if(lv_pm == 0){
        plotting_results(ME_levels, cat, MAR_rel[rel], proportions_mv, napproaches, nconds)
      }else{
        plotting_results(ME_levels, cat, MAR_rel[rel],  proportions_mv, napproaches, nconds, pm_name = lv_pm)
      }
    }
  }
} 

lapply(lv_pm, function(lv_pm){results(ME_levels, lv_pm, MAR_rel, proportions_mv, napproaches, ncat, nconds)})
