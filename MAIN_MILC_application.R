### Description script ### 
# Main script to perform application. 
# First some preparation: load libraries, define parameters
# The main code is divided in 1) apply the different approaches to treat missing data, 
# 2) perform MILC on the data set after missing data was treated
# 3) make tables and figures from the results, 
# 4) perform all previous defined results to obtain the tables and figures

### Load libraries ### 
library(poLCA)
library(parallel)
library(doParallel)
library(doRNG)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mice)

### define parameters ###
nboots = c(10, 5)
ncat = 7
napproaches = 3
data_sets = c("gas","electricity")
pms = c("lv", "bias", "RMSE" )
seeds = c(1, 2)

########################
### Data preparation ###
########################
# Load data
setwd('//cbsp.nl/Productie/Secundair/MPOnderzoek/Werk/Combineren/Projecten/MeetfoutenEnergiebestanden/Data/bestand 3/2020/')
data = readRDS(file = 'data6bron_aangepast.rds')
# Set working directory in own map
setwd('//cbsp.nl/Productie/Secundair/MPOnderzoek/Werk/Combineren/Projecten/MeetfoutenEnergiebestanden/applicatie_iris')
data$energiedrager    = as.factor(data$energiedrager)
data$ABRP8[data$ABRP8 == ""]                       = NA
data$Dataland8[data$Dataland8 == ""]               = NA
data$Kernwoord8[data$Kernwoord8 == ""]             = NA
data$Definitief8[data$Definitief8 == ""]           = NA

data$ABRP8       = data$ABRP8      [, drop = T]
data$Dataland8   = data$Dataland8  [, drop = T]
data$Kernwoord8  = data$Kernwoord8 [, drop = T]
data$Definitief8 = data$Definitief8[, drop = T]

##########################
### Treat missing data ###
##########################
apply_approach = function(data, approach, seed){
  # small function
  convert_na_to_factor <- function(x) {
    y = as.numeric(x)
    m = max(y, na.rm = T)
    y[is.na(y)] = m + 1
    y
  }
  if(seed == 1){
    set.seed(123)
  }else if(seed == 2){
    set.seed(987)
  }
  samp_data = data[sample(0.1*nrow(data), replace =F), ]
  if(approach == 1){# CCA
    data_approach = samp_data[complete.cases(samp_data[c("ABRP8", "Dataland8", "Kernwoord8")]), ]
    for (column in c('ABRP8', 'Dataland8', 'Kernwoord8')) {
      data_approach[, column] = as.numeric(data_approach[, column])
    }
  }else if(approach == 2){# ECMV
    data_approach = samp_data
    for(column in c('ABRP8', 'Dataland8', 'Kernwoord8')){
      data_approach[, column] = convert_na_to_factor(data_approach[, column])
    }
  }else if(approach == 3){ # MICE
    ind = samp_data[,c(12, 13, 16)]
    log = capture.output({ MICE = mice(ind, method = "polyreg")})
    imp = ind
    imp = complete(MICE)
    for (column in c('ABRP8', 'Dataland8', 'Kernwoord8')) {
      imp[, column] = as.numeric(imp[, column])
    }
    data_approach = samp_data
    data_approach$ABRP8 = imp$ABRP8
    data_approach$Dataland8 = imp$Dataland8
    data_approach$Kernwoord8 = imp$Kernwoord8
  }
  #Use this data set
  cleaned_data =  data_approach
  y = cleaned_data[,c(12, 13, 16)] 
  n = nrow(cleaned_data)
  return(list("cleaned_data" = cleaned_data, "y" = y, "n" = n))
}

##########################
### Apply MILC ###
##########################
apply_MILC = function(cleaned_data, data_type, y, n, nboot, ncat, seed){
  res = matrix(data = NA, nrow = n, ncol = nboot)  # storing imputed values
  res = as.data.frame(res)
  names(res) = c(1:nboot)
  theta_LC =  matrix(data = NA, nrow = nboot, ncol = ncat) # storing each imputation's parameter
  var_LC = matrix(data = NA, nrow = nboot, ncol = ncat) # storing each imputation's parameter
  # Start bootstrap and perform MILC
  starttime=Sys.time()
  for(boot in 1:nboot){
    cat("boot =",boot,"\n\n")
    set.seed(boot)
    Ready = T
    while(Ready){
      Bsample = cleaned_data[sample(1:n, n, replace = T), ]
      
      ready=FALSE
      while(!ready){
        LCA = tryCatch(expr = poLCA(formula = cbind(ABRP8, Dataland8, Kernwoord8) ~ 1, 
                                    data = Bsample, nclass = ncat, verbose = T, maxiter = 10000,
                                    nrep = 10, calc.se =F), error = function(e) NULL)
        ready = !is.null(LCA)
      }
      
      K = NA
      for(k in 1:ncat){
        K[k] = which.max(LCA$probs$ABRP8[, k])
      }
      
      K = as.factor(K)
      prob.start.new = poLCA.reorder(LCA$probs.start, K)
      
      LCAr = tryCatch(expr = poLCA(formula = cbind(ABRP8, Dataland8, Kernwoord8) ~ 1, 
                                   data = Bsample, nclass = ncat, verbose = T, 
                                   probs.start = prob.start.new, maxiter = 10000, calc.se =F),
                      error = function(e) NULL)
      probs = LCAr$probs
      write.csv2(probs, file = sprintf("probs_%s_approach%s_boot%s_%sb_seed%s.csv", data_type, approach, boot, nboot, seed), quote =F, row.names = F)
      ppd = tryCatch(expr = poLCA.posterior(LCAr, y = y), error = function(e) NA)
      res[ , boot] = tryCatch(expr = rmulti(ppd), error = function(e) NA)
      for (cat in 1:ncat){
        theta_LC[boot,cat] = sum(res[,boot] == cat, na.rm =T)/nrow(res)
        var_LC[boot,cat] = theta_LC[boot,cat]*(1-theta_LC[boot,cat])/nrow(res) 
      }
      Ready = any(is.na(var_LC[boot,]))
      if(Ready){
        errorrecord = c(errorrecord,boot)
      }
    }
  }
  SimulationTotalTime = Sys.time()-starttime
  print(SimulationTotalTime)
  return(list("theta_LC" = theta_LC, "var_LC" = var_LC))
}

###############################
### Obtain and plot results ###
###############################
output_MILC = function(theta_LC, var_LC, theta, data_type, approach, nboot, ncat, seed){
  ### Obtain  results ###
  # 1) LV estimation
  p_mean = apply(theta_LC, MARGIN = 2, FUN = mean) # pooled mean 
  write.csv2(p_mean, file = sprintf("mean_%s_approach%s_%sb_seed%s.csv", data_type, approach, nboot, seed), quote =F, row.names = F)
  # 2) Bias estimation
  bias = data.frame(t(theta_LC))
  theta_df = data.frame(theta)
  for(col in 1:ncol(bias)){
    bias[ , col] = bias[,col] - theta_df
  }
  p_bias = apply(bias, MARGIN = 1, FUN = mean) # pooled bias
  write.csv2(p_bias, file = sprintf("bias_%s_approach%s_%sb_seed%s.csv", data_type, approach, nboot, seed), quote =F, row.names = F)
  # 3) Variance
  theta_LC_df = data.frame(t(theta_LC))
  V_within = apply( var_LC/nboot, MARGIN = 2, FUN = mean) 
  V_between = apply((theta_LC_df-p_mean)^2, MARGIN = 1, FUN = sum)/(nboot-1) 
  p_VAR = V_within + V_between + (V_between/nboot)
  # 4) MSE
  p_MSE = p_VAR + p_bias^2
  # 5) RMSE
  p_RMSE = sqrt(p_MSE)
  write.csv2(p_RMSE, file = sprintf("RMSE_%s_approach%s_%sb_seed%s.csv", data_type, approach, nboot, seed), quote=FALSE, row.names = F)
  ### To plot results ###
  plot.data = data.frame(category = rep(1:ncat),
                         lv = p_mean,
                         bias = p_bias,
                         RMSE = p_RMSE,
                         sd = sqrt(p_VAR)*1.96)
  write.csv2(plot.data, file = sprintf("plot_table_%s_approach%s_%sb_seed%s.csv", data_type, approach, nboot, seed), quote =F, row.names = F)
}

plot_MILC = function(pms, data_type, approach, theta, nboot, ncat, seed){
  new_x = c('A', 'BCD', 'EF', 'GHI', 'JKLMN', 'OPQ', 'RSTU')
  for(pm in 1:length(pms)){
    if(pm == 2){
      if(data_type == "electricity"){
        y_lim = c(-0.2, 0.4)
      }else if(data_type == "gas"){
        y_lim = c(-0.25, 0.6)
      }
    }else if(pm == 3){
      if(data_type == "electricity"){
        y_lim = c(0, 0.4)
      }else if(data_type == "gas"){
        y_lim = c(0, 0.45)
      }
    }else if(pm == 1){
      if(data_type == "electricity"){
        y_lim = c(0, 0.6)
      }else if(data_type == "gas"){
        y_lim = c(-0.2, 0.6)
      }
    }
    table_pm = read.csv(file=sprintf("plot_table_%s_approach%s_%sb_seed%s.csv", data_type, approach, nboot, seed), header = T, sep = ";")
    table_pm[1:5] = lapply(table_pm[1:5], gsub, pattern = ",", replacement = ".")
    table_pm[] = as.numeric(as.matrix(table_pm))
    if(pm == 1){
      ggplot(data = table_pm, aes(x = category, y = table_pm[, 1+pm], ymax = lv+sd, ymin = lv-sd, colour = "red")) +
        geom_point(aes(x =  1:ncat, y = theta), col = 'black', shape = 17, size = 1.5)+
        geom_errorbar(width = 0.2, position = position_dodge(0.03)) +
        geom_point() +
        theme_minimal() +
        theme(legend.position ="none") +
        ylim(y_lim) +
        scale_x_continuous(breaks = c(1:ncat), labels = new_x ) +
        ylab("Estimate") +
        xlab("")
      ggsave(sprintf("plot_%s_%s_approach%s_%sb_seed%s.pdf", pms[pm], data_type, approach, nboot, seed), width = 5, height = 5)
    }else{
      pms = c("lv", "Bias", "RMSE" )
      ggplot(data = table_pm, aes(x = category, y = table_pm[, 1+pm], colour = "red")) +
        geom_point() +
        theme_minimal() +
        theme(legend.position ="none") +
        theme(text = element_text(size = 15)) +
        #ylim(y_lim) +
        scale_x_continuous(breaks = c(1:ncat), labels = new_x ) +
        ylab(sprintf("%s", pms[pm]) ) +
        xlab("")
      pms = c("lv", "bias", "RMSE" )
      ggsave(sprintf("plot_%s_%s_approach%s_%sb_seed%s.pdf", pms[pm], data_type, approach, nboot, seed),width = 5, height = 5)
   }
  }
}

#############################
### Perm simulation study ###
#############################
for(seed in seeds){
  print(seed)
  for(nboot in nboots){
    print(nboot)
    for(dat in data_sets){
      if(dat == "electricity"){
        data_set = data[data$energiedrager == "E",]
      }else if(dat == "gas"){
        data_set = data[data$energiedrager == "G",]
      }
      perc_mv =(colMeans(is.na(data_set)))*100
      write.csv2(perc_mv, file=sprintf("perc_mv_%s_%sb_%s.csv", dat, nboot, seed), quote=FALSE, row.names=F)
      n_data = nrow(data_set)
      write.csv2(n_data, file=sprintf("n_data_%s_%sb_%s.csv", dat, nboot, seed), quote=FALSE, row.names=F)
      theta = matrix(data = NA, nrow = ncat, ncol = 1)
      for(cat in 1:ncat){
        theta[cat] = table(data_set$Definitief8)[cat] 
      }
      theta = theta/sum(theta)
      write.csv2(theta, file=sprintf("theta_%s_%sb_%s.csv", dat, nboot, seed), quote=FALSE, row.names=F)
      for(approach in 1:napproaches){
        print(approach)
        results_approach = apply_approach(data_set, approach, seed)
        write.csv2(results_approach$n, file=sprintf("n_approach%s_%s_%sb_%s.csv", approach, dat, nboot, seed), quote =F, row.names=F)
        results_MILC = apply_MILC(results_approach$cleaned_data, dat, results_approach$y, results_approach$n, nboot, ncat, seed)
        output = output_MILC(results_MILC$theta_LC, var_LC = results_MILC$var_LC, theta, dat, approach, nboot, ncat, seed)
        result = plot_MILC(pms, dat, approach, theta, nboot, ncat, seed)
      }
    }
  }
}


# If you only want to adjust the plots
# for(seed in seeds){
#   print(seed)
#   for(nboot in nboots){
#     print(nboot)
#     for(dat in data_sets){
#       if(dat == "electricity"){
#         data_set = data[data$energiedrager == "E",]
#       }else if(dat == "gas"){
#         data_set = data[data$energiedrager == "G",]
#       }
#       
#       theta = matrix(data = NA, nrow = ncat, ncol = 1)
#       for(cat in 1:ncat){
#         theta[cat] = table(data_set$Definitief8)[cat] 
#       }
#       theta = theta/sum(theta)
#       for(approach in 1:napproaches){
#         print(approach)
#         result = plot_MILC(pms, dat, approach, theta, nboot, ncat, seed)
#       }
#     }
#   }
# }



