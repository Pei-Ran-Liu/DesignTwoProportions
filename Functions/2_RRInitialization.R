################################################################
##### Initialization of the Rejection Region (Algorithm 1) #####
################################################################
library(foreach)
library(doParallel)

##### For Frequentist Methods #####
binary_search_freq_y1 <- function(y2, n1, n2, delta, gamma, method="Unpooled Z"){
  # y2: for a given y2, do binary search of y1
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # gamma: cut value of the decision -  statistic > gamma
  # method: which statistic used
  
  # check the method
  if(method=="Pooled Z"){
    func <- z_pool
  } else if(method=="Unpooled Z"){
    func <- z_unpool
  } else if(method=="score"){
    func <- score
  } else if(method=="Fisher's condition"){
    func <- fisher_condp
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  L <- 0
  U <- n1
  t_L <- func(L,y2,n1,n2,delta)
  t_U <- func(U,y2,n1,n2,delta)
  if(t_L>gamma){
    return(L)
  } else if(t_U<=gamma){
    return(n1+1)
  } else {
    while(U-L>1){
      M <- ceiling((L+U)/2)
      t_M <- func(M,y2,n1,n2,delta)
      if(t_M>gamma){
        U <- M
      } else{
        L <- M
      }
    }
    return(U)
  }
  
}

binary_search_freq_y2 <- function(y1, n1, n2, delta, gamma, method="Unpooled Z"){
  # y1: for a given y1, do binary search of y2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # gamma: cut value of the decision -  statistic > gamma
  # method: which statistic used
  
  # check the method
  if(method=="Pooled Z"){
    func <- z_pool
  } else if(method=="Unpooled Z"){
    func <- z_unpool
  } else if(method=="score"){
    func <- score
  } else if(method=="Fisher's condition"){
    func <- fisher_condp
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  L <- 0
  U <- n2
  t_L <- func(y1,L,n1,n2,delta) 
  t_U <- func(y1,U,n1,n2,delta) 
  if(t_L<=gamma){
    return(-1)
  } else if(t_U>gamma){
    return(n2)
  } else {
    while(U-L>1){
      M <- floor((L+U)/2)
      t_M <- func(y1,M,n1,n2,delta)
      if(t_M>gamma){
        L <- M
      } else{
        U <- M
      }
    }
    return(L)
  }
  
}

# combine all together 
freq_rr <- function(n1, n2, delta, gamma, method="Unpooled Z"){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # gamma: cut value of the decision -  statistic > gamma
  # method: which statistic used
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  if( ((delta>0)&(method=="Unpooled Z"))||(n2>n1) ){
    return(foreach (i=seq(0,n1,1), .combine=c) %dopar% {
      binary_search_freq_y2(i, n1, n2, delta, gamma, method)
    })
  } else{
    return(foreach (i=seq(0,n2,1), .combine=c) %dopar% {
      binary_search_freq_y1(i, n1, n2, delta, gamma, method)
    })
  }
  
  stopImplicitCluster()
  
}



##### For Bayesian Probability Methods #####
binary_search_bayes_y1 <- function(y2, n1, n2, delta, gamma, method="Independent Beta", parameters=c(1,1,1,1)){
  # y2: for a given y2, do binary search of y1
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # gamma: cut value of the decision -  statistic > gamma
  # method = "", and parameters of the corresponding method = parameters
  # Independent Beta (two independent beta prior with parameters a1, b1, a2, b2)
  # Logit Normal (logit-normal prior with parameters mu1, sigma1 (sd), mu2, sigma2(sd), sigma12 (cov))
  
  # check the method
  if(method=="Independent Beta"){
    func <- bayes_indp_beta
  } else if(method=="Logit Normal"){
    func <- bayes_logit
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  L <- 0
  U <- n1
  t_L <- func(L,y2,n1,n2,delta,parameters)
  t_U <- func(U,y2,n1,n2,delta,parameters)
  if(t_L>gamma){
    return(L)
  } else if(t_U<=gamma){
    return(n1+1)
  } else {
    while(U-L>1){
      M <- ceiling((L+U)/2)
      t_M <- func(M,y2,n1,n2,delta,parameters)
      if(t_M>gamma){
        U <- M
      } else{
        L <- M
      }
    }
    return(U)
  }
  
}

binary_search_bayes_y2 <- function(y1, n1, n2, delta, gamma, method="Independent Beta", parameters=c(1,1,1,1)){
  # y1: for a given y1, do binary search of y2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # gamma: cut value of the decision -  statistic > gamma
  # method = "", and parameters of the corresponding method = parameters
  # Independent Beta (two independent beta prior with parameters a1, b1, a2, b2)
  # Logit Normal (logit-normal prior with parameters mu1, sigma1 (sd), mu2, sigma2(sd), sigma12 (cov))
  
  # check the method
  if(method=="Independent Beta"){
    func <- bayes_indp_beta
  } else if(method=="Logit Normal"){
    func <- bayes_logit
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  L <- 0
  U <- n2
  t_L <- func(y1,L,n1,n2,delta,parameters)
  t_U <- func(y1,U,n1,n2,delta,parameters)
  if(t_L<=gamma){
    return(-1)
  } else if(t_U>gamma){
    return(n2)
  } else {
    while(U-L>1){
      M <- floor((L+U)/2)
      t_M <- func(y1,M,n1,n2,delta,parameters)
      if(t_M>gamma){
        L <- M
      } else{
        U <- M
      }
    }
    return(L)
  }
  
}

# combine all together
bayes_rr <- function(n1, n2, delta, gamma, method="Independent Beta", parameters=c(1,1,1,1)){
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  if(n1>=n2){
    return(foreach (i=seq(0,n2,1), .combine=c) %dopar% {
      binary_search_bayes_y1(i, n1, n2, delta, gamma, method, parameters)
    })
  } else{
    return(foreach (i=seq(0,n1,1), .combine=c) %dopar% {
      binary_search_bayes_y2(i, n1, n2, delta, gamma, method, parameters)
    })
  }
  
  stopImplicitCluster()
  
}




