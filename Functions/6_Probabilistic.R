##############################################################
##### Randomized Part of the Corresponding Decision Rule #####
##############################################################

# Step 1: generate the second part of the critical region corresponding to T(y_1,y_2)==gamma
secondp_rd_freq <- function(boundary_deter, n1, n2, delta, method="Unpooled Z"){
  
  if(((length(boundary_deter)==n1+1)&(n1!=n2))){
    
    candidate <- boundary_deter + 1
    tsvalue_candidate <- freq_ts_value(candidate, n1, n2, delta, method)
    tsvalue_candidate[candidate>n2] <- -Inf
    gamma <-  max(tsvalue_candidate)
    res_select <- which(tsvalue_candidate==gamma)
    
  } else{
    
    candidate <- boundary_deter - 1
    tsvalue_candidate <- freq_ts_value(candidate, n1, n2, delta, method)
    tsvalue_candidate[candidate<0] <- -Inf
    gamma <-  max(tsvalue_candidate)
    res_select <- which(tsvalue_candidate==gamma)
    
  }
  
  return(res_select)
  
}

secondp_rd_bayes <- function(boundary_deter, n1, n2, delta, method="Independent Beta", parameters=c(1,1,1,1)){
  
  if(((length(boundary_deter)==n1+1)&(n1!=n2))){
    
    candidate <- boundary_deter + 1
    tsvalue_candidate <- bayes_ts_value(candidate, n1, n2, delta, method, parameters)
    tsvalue_candidate[candidate>n2] <- -Inf
    gamma <-  max(tsvalue_candidate)
    res_select <- which(tsvalue_candidate==gamma)
    
  } else{
    
    candidate <- boundary_deter - 1
    tsvalue_candidate <- bayes_ts_value(candidate, n1, n2, delta, method, parameters)
    tsvalue_candidate[candidate<0] <- -Inf
    gamma <-  max(tsvalue_candidate)
    res_select <- which(tsvalue_candidate==gamma)
    
  }
  
  return(res_select)
  
}

# Step 2. generate the expected value of the part obtained at step 1
powerf_2nd <- function(boundary_deter, res_select, n1, n2, delta, theta){
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  if(((length(boundary_deter)==n1+1)&(n1!=n2))){
    
    candidate <- boundary_deter + 1
    return(foreach(i=theta, .combine=c) %dopar% {
      sum(exp(sapply(res_select-1, dbinom, size=n1, prob=i, log = TRUE) + sapply(candidate[res_select], dbinom, size=n2, prob=i+delta, log = TRUE)))
    })
    
  } else{
    
    candidate <- boundary_deter - 1
    return(foreach(i=theta, .combine=c) %dopar% {
      sum(exp(sapply(res_select-1, dbinom, size=n2, prob=i+delta, log = TRUE) + sapply(candidate[res_select], dbinom, size=n1, prob=i, log = TRUE)))
    })
    
  }
  
  stopImplicitCluster()
  
}

# Step 3. calculate the proper value of xi such that the actual type I error is equal to the intended one exactly
binarys_xi <- function(part1_tIe, part2_tIe, alpha, precision = 0.00001){
  xi_L <- 0
  tIe_L <- max(part1_tIe)
  xi_U <- 1
  tIe_U <- max(part1_tIe + part2_tIe)
  
  if(max(part1_tIe)>alpha){
    
    print("The deterministic part is not specified in a right way!")
    return(NA)
    
  } else{
    
    while( ((xi_U-xi_L)>precision || (round(tIe_L, -4*log10(precision))<alpha)) ){ # -log10(precision) # 1/precision
      xi_M <- (xi_L+xi_U)/2
      tIe_M <- max(part1_tIe + xi_M*part2_tIe)
      if(tIe_M <= alpha){
        xi_L <- xi_M
        tIe_L <- tIe_M
      } else{
        xi_U <- xi_M
        tIe_U <- tIe_M
      }
    }
    return(c(xi_L, tIe_L))
    
  }
  
}


