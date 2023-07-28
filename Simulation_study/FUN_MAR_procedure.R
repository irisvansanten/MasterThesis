### Description script ###
# Example script on how the scaled probabilities and Cramérs V were calculated when MAR mechanisms is employed. 
# Due to many simulations this is an approximation of the values in the simulation study.
# First observed variable Z will be simulated and then scaled probabilities are calculated 

### Load libraries ###
library(Rlab)
library(confintr)
library(poLCA)

### Define parameters ###
k = 2 
nsize = 10000
ncat = 4
cat_miss_prop_strong = c(0.2, 0.5, 0.01, 0.6)
cat_miss_prop_weak = c(0.4, 0.3, 0.35, 0.3)
T_dist = c(0.1, 0.2, 0.3, 0.4)
total_miss_prop = c(0.1, 0.2, 0.3, 0.4)

##################
### Simulate Z ###
##################

set.seed(999)
model = poLCA.simdata(N       = nsize,
                      nclass  = ncat,
                      probs   = list(matrix(c(285/300, 5/300, 5/300, 5/300,
                                              5/300, 285/300, 5/300, 5/300, 
                                              5/300, 5/300, 285/300, 5/300,
                                              5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow = T),
                                     matrix(c(285/300, 5/300, 5/300, 5/300,
                                              5/300, 285/300, 5/300, 5/300, 
                                              5/300, 5/300, 285/300, 5/300,
                                              5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow = T),
                                     matrix(c(285/300, 5/300, 5/300, 5/300,
                                              5/300, 285/300, 5/300, 5/300, 
                                              5/300, 5/300, 285/300, 5/300,
                                              5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow = T),
                                     matrix(c(285/300, 5/300, 5/300, 5/300,
                                              5/300, 285/300, 5/300, 5/300, 
                                              5/300, 5/300, 285/300, 5/300,
                                              5/300, 5/300, 5/300, 285/300), ncol = ncat, byrow = T),
                                     matrix(c(18/30, 4/30, 4/30, 4/30,
                                              8/30, 6/30, 8/30, 8/30, 
                                              2/30, 2/30, 24/30, 2/30,
                                              5/30, 5/30, 5/30, 15/30), ncol = ncat, byrow = T) ),
                      P       = T_dist,
                      missval = F)
df_Z = as.data.frame(model$dat$Y5)
known_prop = table(unlist(df_Z))/(sum(table(unlist(df_Z))))
print("known_prop")
print(known_prop)

#############################################################
### Calculating Scaled probabilities and Cramér's V value ###
#############################################################

MAR_procedure = function(known_prop, total_miss_prop, cat_miss_prop, nsize, k){
  weight = total_miss_prop/(sum(known_prop*cat_miss_prop))
  cat_miss = weight*known_prop*cat_miss_prop
  cat_obs = known_prop - cat_miss
  sc = cat_miss/known_prop
  print("scaled probability")
  print(sc)
  # To calculate relation 
  keeps = c("Freq", "Freq.1")
  df = data.frame(cat_miss, cat_obs)
  contigency_table = df[keeps]*nsize
  print("Contigency table")
  print(contigency_table)
  test = chisq.test(contigency_table)
  N_strong = sum(contigency_table)
  cramersv = sqrt(test$statistic/(N_strong*(k-1)))
  print("Cramér's V value")
  print(cramersv)
}

# For strong MAR mechanisms
for(miss_prop in total_miss_prop){
  MAR_procedure(known_prop, miss_prop, cat_miss_prop_strong, nsize, k)
}

# For weak MAR mechanisms
for(miss_prop in total_miss_prop){
  MAR_procedure(known_prop, miss_prop, cat_miss_prop_weak, nsize, k)
}