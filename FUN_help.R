### Description script ###
# Functions that are used when missing values are introduced and treated, and when MILC was performed.
# These functions are used in the code when we parallize. 

#############################
### Missing data function ###
#############################

# 1) Employ MAR mechanism
MAR_function = function(df, cat_miss_prop, proportions_mv, ncat){
  df_NA = data.frame(matrix(nrow = nrow(df), ncol = ncol(df)))
  colnames(df_NA) = c("Y1", "Y2", "Y3", "Y4", "Y5" )
  known_prop = table(unlist(df[, 1:4]))/(sum(table(unlist(df[, 1:4]))))
  weight = proportions_mv/(sum(known_prop*cat_miss_prop))
  cat_miss = weight*known_prop*cat_miss_prop
  cat_miss_prob = cat_miss / known_prop
  for (row in 1:nrow(df)){
    row_mv = df[row,]
    if (row_mv$Y5 == 1 ){
      random_values = append(rbern(ncat, prob = 1-cat_miss_prob[row_mv$Y5]), row_mv$Y5)
      df_NA[row,] = random_values
    }else if(row_mv$Y5 == 2){
      random_values = append(rbern(ncat, prob = 1-cat_miss_prob[row_mv$Y5]), row_mv$Y5)
      df_NA[row,] = random_values
    }else if(row_mv$Y5 == 3){
      random_values = append(rbern(ncat, prob = 1-cat_miss_prob[row_mv$Y5]), row_mv$Y5)
      df_NA[row,] = random_values
    }else if(row_mv$Y5 == 4){
      random_values = append(rbern(ncat, prob = 1-cat_miss_prob[row_mv$Y5]), row_mv$Y5)
      df_NA[row,] = random_values
    }
    df_NA[df_NA == 0] = NA
    df_MAR = df
    df_MAR[which(is.na(df_NA), arr.ind = T)] = NA
  }
  return(df_MAR)
}

# 2) Introduce missing values with MAR mechanism
missing_values = function(model, 
                          cat_miss_prop, 
                          proportions_mv,
                          ncat
){
  dataframe_model_NA = MAR_function(model$dat, cat_miss_prop, proportions_mv, ncat)
  dataframe_NA = cbind(dataframe_model_NA, Tclass = model$Tclass)
  return(list("dataframe_NA" = dataframe_NA, "dataframe_model_NA" = dataframe_model_NA, "model" = model)) 
}

# 3) Treat missing values with 3 different approaches
# Approach 1: Omit the rows with missing values
# Approach 2: Assign the missing values to an extra category
# Approach 3: Impute the missing values with mode and MICE imputation
approaches_function= function(dataframe_NA, dataframe_model_NA, model, approach){
  if(approach == 1){
    dataframe_mv= na.omit(dataframe_NA)
  }else if(approach == 2){
    dataframe_NA[is.na(dataframe_NA)] = as.factor(5)
    dataframe_mv = dataframe_NA
  }else if(approach == 3){
    log = capture.output({ MICE = mice(dataframe_model_NA, method = "polyreg")})
    dataframe_NA_imp = dataframe_model_NA
    dataframe_NA_imp = complete(MICE)
    dataframe_mv = cbind(dataframe_NA_imp, Tclass = model$Tclass)
  }
  return(dataframe_mv)
}

#####################
### MILC function ###
#####################

# 4) Bootstrap data
bootdata = function(data, 
                    nboot){
  boots = rmultinom(nboot, 
                    sum(data$n), 
                    data$n/sum(data$n))
  dfboot = cbind(data, boots)
  names(dfboot) = c(names(data), paste0('b', 1:nboot))
  return(dfboot)
}


# 5) Posteriors function in which LC analysis is performed and latent class is imputed. 
# LC model consists of four indicators (Y1 to Y4) with the same amount of measurement error. 
posteriors_no_Z = function(bootdat, k){
  log = capture.output({ 
    longdat = as.data.frame(fre2dat(bootdat[, c(1:5,k+6)]))                     
  })
  LCA =  poLCA(formula = cbind(Y1, Y2, Y3, Y4) ~ 1,  
                data    = longdat,
                nclass  = 4,      
                nrep    = 10)
  order  =  c(which.max(LCA$probs$Y1[, 1]),                                     
               which.max(LCA$probs$Y1[, 2]), 
               which.max(LCA$probs$Y1[, 3]), 
               which.max(LCA$probs$Y1[, 4]))
  probs.start.new = poLCA.reorder(LCA$probs.start, order)                      
  LCAr    = poLCA(formula = cbind(Y1, Y2, Y3, Y4) ~ 1,
                   data    = longdat,
                   nclass  = 4,
                   probs.start = probs.start.new)
  entr   = poLCA.entropy(LCAr)  
  ordat  =  as.data.frame(fre2dat(bootdat[, c(1:6)]))
  posdat =  cbind(ordat, 
                   poLCA.posterior(lc = LCAr, y = ordat[, -5]),                 
                   entr = entr) 
  return(posdat)
}

# LC model consists of five indicators for which Y1 to Y4 have the same amount of measurement errors 
# and Y5 a different and higher amount of measurement errors.  
posteriors_Z = function(bootdat, k){
  log = capture.output({ 
    longdat = as.data.frame(fre2dat(bootdat[,c(1:6, k+7)]))                     
  })
  LCA =  poLCA(formula = cbind(Y1, Y2, Y3, Y4, Y5) ~ 1,  
                data    = longdat,
                nclass  = 4,      
                nrep    = 10)
  order  =  c(which.max(LCA$probs$Y1[, 1]),                                     
               which.max(LCA$probs$Y1[, 2]), 
               which.max(LCA$probs$Y1[, 3]), 
               which.max(LCA$probs$Y1[, 4]))
  probs.start.new = poLCA.reorder(LCA$probs.start, order)                      
  LCAr    = poLCA(formula = cbind(Y1, Y2, Y3, Y4, Y5) ~ 1,
                   data    = longdat,
                   nclass  = 4,
                   probs.start = probs.start.new)
  entr   = poLCA.entropy(LCAr)  
  ordat  =  as.data.frame(fre2dat(bootdat[, c(1:7)]))
  
  posdat =  cbind(ordat, 
                   poLCA.posterior(lc = LCAr, y = ordat[, -6]),                 
                   entr = entr) 
  return(posdat)
}