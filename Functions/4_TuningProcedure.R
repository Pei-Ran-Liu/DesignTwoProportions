################################################################
##### Tuning procedure (Forward and Backward Algorithms 2) #####
################################################################

##### Forward Procedure (Freq) #####
forward_freq <- function(n1, n2, delta, alpha, 
                         boundary_initial, type_I_error_initial, 
                         theta, method="Unpooled Z"){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 -delta
  #  -1 < delta < 1
  # alpha: intended significance level
  # boundary_initial (vector): boundary of the rejection region obtained based on the initial gamma
  # type_I_error_initial (vector): type I error at each grid based on the initial gamma
  # theta (vector): boundary of p1 with p2=p1+delta
  
  indicator <- 1
  
  # select the test statistic
  if(method=="Pooled Z"){
    func <- z_pool
  } else if(method=="Unpooled Z"){
    func <- z_unpool
  } else if(method=="score"){
    func <- score
  } else if(method=="Fisher's condition"){
    func <- fisher_condp
  } else{
    print("This method is not available now!")
    return(NA)
  }
  
  if( (((delta>0)&(method=="Unpooled Z")) || (n2>n1)) ){
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    boundary_value <- freq_ts_value(boundary, n1, n2, delta, method)
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    while((max(type_I_error)>alpha)&indicator){
      ind_null <- which.max(type_I_error)
      while(type_I_error[ind_null]>alpha){
        gamma <- min(boundary_value)
        ind <- which(boundary_value==gamma) 
        n_ind <- length(ind)
        for(i in 1:n_ind){
          type_I_error <- type_I_error - exp(sapply(theta+delta, dbinom, x=boundary[ind[i]], size=n2, log=T)) * 
            exp(sapply(theta, dbinom, x=ind[i]-1, size=n1, log=T))
          boundary[ind[i]] <- boundary[ind[i]]-1
          if(boundary[ind[i]]<0){
            boundary_value[ind[i]] <- Inf
            if(min(boundary_value)==Inf){
              indicator <- 0
              break;
            }
          } else {
            boundary_value[ind[i]] <- func(ind[i]-1,boundary[ind[i]],n1,n2,delta)
          }
        }
        iteration <- iteration + 1
      }
    }
    if(indicator==0){
      print("Rejection region is the domain (reject the null all the time)")
      boundary <- rep(-1,n1+1)
      return(list(gamma, boundary, 0, iteration, boundary_initial, type_I_error_initial))
    } else{
      return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
    }
    
  } else{ 
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    boundary_value <- freq_ts_value(boundary, n1, n2, delta, method)
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    while((max(type_I_error)>alpha)&indicator){
      ind_null <- which.max(type_I_error)
      while(type_I_error[ind_null]>alpha){
        gamma <- min(boundary_value)
        ind <- which(boundary_value==gamma) 
        n_ind <- length(ind)
        for(i in 1:n_ind){
          type_I_error <- type_I_error - exp(sapply(theta, dbinom, x=boundary[ind[i]], size=n1, log=T)) * 
            exp(sapply(theta+delta, dbinom, x=ind[i]-1, size=n2, log=T))
          boundary[ind[i]] <- boundary[ind[i]]+1
          if(boundary[ind[i]]>n1){
            boundary_value[ind[i]] <- Inf
            if(min(boundary_value)==Inf){
              indicator <- 0
              break
            }
          } else {
            boundary_value[ind[i]] <- func(boundary[ind[i]],ind[i]-1,n1,n2,delta)
          }
        }
        iteration <- iteration + 1
      }
    }
    if(indicator==0){
      print("Rejection region is the domain (reject the null all the time)")
      boundary <- rep(n1+1,n2+1)
      return(list(gamma, boundary, 0, iteration, boundary_initial, type_I_error_initial))
    } else{
      return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
    }
    
  }
  
}



##### Backward Procedure (Freq) #####
backward_freq <- function(n1, n2, delta, alpha, 
                          boundary_initial, type_I_error_initial, 
                          theta, method="Unpooled Z"){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin with null: p1 <= p2 -delta
  #  -1 < delta < 1
  # alpha: intended significance level
  # boundary_initial (vector): boundary of the rejection region obtained based on the initial gamma
  # type_I_error_initial (vector): type I error at each grid based on the initial gamma
  # theta (vector): boundary of p1 with p2=p1+delta
  
  # select the test statistic
  if(method=="Pooled Z"){
    func <- z_pool
  } else if(method=="Unpooled Z"){
    func <- z_unpool
  } else if(method=="score"){
    func <- score
  } else if(method=="Fisher's condition"){
    func <- fisher_condp
  } else{
    print("This method is not available now!")
    return(NA)
  }
  
  if( (((delta>0)&(method=="Unpooled Z")) || (n2>n1)) ){
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    candidate <- boundary_initial + 1
    boundary_value <- freq_ts_value(candidate, n1, n2, delta, method)
    boundary_value[candidate>n2] <- -Inf
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    indicator <- 1 
    while(indicator==1){
      gamma_temp <- max(boundary_value)
      ind_temp <- which(boundary_value==gamma_temp)
      n_temp <- length(ind_temp)
      type_I_error_temp <- type_I_error
      
      for(i in 1:n_temp){
        type_I_error_temp <- type_I_error_temp + exp(sapply(theta+delta, dbinom, x=candidate[ind_temp[i]], size=n2, log=T)) * 
          exp(sapply(theta, dbinom, x=ind_temp[i]-1, size=n1, log=T)) 
        candidate[ind_temp[i]] <- candidate[ind_temp[i]]+1
        if(candidate[ind_temp[i]]>n2){
          boundary_value[ind_temp[i]] <- -Inf
        } else {
          boundary_value[ind_temp[i]] <- func(ind_temp[i]-1,candidate[ind_temp[i]],n1,n2,delta)
        }
      }
      
      if(max(type_I_error_temp)>alpha){
        indicator <- 0 
        if(iteration==0){
          print("gamma does not change")
          return(list("No", boundary_initial, type_I_error_initial, 0, boundary_initial, type_I_error_initial))
        } else{
          gamma <- gamma_temp
          return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
        }
      } else{
        boundary[ind_temp] <- boundary[ind_temp] + 1
        type_I_error <- type_I_error_temp
      }
      iteration <- iteration + 1
    }
    
  } else{
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    candidate <- boundary_initial - 1
    boundary_value <- freq_ts_value(candidate, n1, n2, delta, method)
    boundary_value[candidate<0] <- -Inf
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    indicator <- 1
    while(indicator==1){
      gamma_temp <- max(boundary_value)
      ind_temp <- which(boundary_value==gamma_temp)
      n_temp <- length(ind_temp)
      type_I_error_temp <- type_I_error
      
      for(i in 1:n_temp){
        type_I_error_temp <- type_I_error_temp + exp(sapply(theta, dbinom, x=candidate[ind_temp[i]], size=n1, log=T)) * 
          exp(sapply(theta+delta, dbinom, x=ind_temp[i]-1, size=n2, log=T)) 
        candidate[ind_temp[i]] <- candidate[ind_temp[i]]-1
        if(candidate[ind_temp[i]]<0){
          boundary_value[ind_temp[i]] <- -Inf
        } else {
          boundary_value[ind_temp[i]] <- func(candidate[ind_temp[i]],ind_temp[i]-1,n1,n2,delta)
        }
      }
      
      if(max(type_I_error_temp)>alpha){
        indicator <- 0 
        if(iteration==0){
          print("gamma does not change")
          return(list("No", boundary_initial, type_I_error_initial, 0, boundary_initial, type_I_error_initial))
        } else{
          gamma <- gamma_temp
          return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
        }
      } else{
        boundary[ind_temp] <- boundary[ind_temp] - 1
        type_I_error <- type_I_error_temp
      }
      iteration <- iteration + 1
    }
    
  }
  
}



##### Forward Procedure (Posterior probability) #####
forward_bayes <- function(n1, n2, delta, alpha, 
                          boundary_initial, type_I_error_initial, theta, 
                          method="Independent Beta", parameters=c(1,1,1,1)){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 -delta
  #  -1 < delta < 1
  # alpha: intended significance level
  # boundary_initial (vector): boundary of the rejection region obtained based on the initial gamma
  # type_I_error_initial (vector): type I error at each grid based on the initial gamma
  # theta (vector): boundary of p1 with p2=p1+delta
  # method = "", and parameters of the corresponding method = parameters
  # Independent Beta (two independent beta prior with parameters a1, b1, a2, b2)
  # Logit Normal (logit-normal prior with parameters mu1, sigma1 (sd), mu2, sigma2(sd), sigma12 (cov))
  
  indicator <- 1 
  
  # select the test statistic
  if(method=="Independent Beta"){
    func <- bayes_indp_beta
  } else if(method=="Logit Normal"){
    func <- bayes_logit
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  # select which to update b_{y1}, or a_{y2}
  if( n2>n1 ){
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    boundary_value <- bayes_ts_value(boundary, n1, n2, delta, method, parameters)
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    while((max(type_I_error)>alpha)&indicator){
      ind_null <- which.max(type_I_error)
      while(type_I_error[ind_null]>alpha){
        gamma <- min(boundary_value)
        ind <- which(boundary_value==gamma) 
        n_ind <- length(ind)
        for(i in 1:n_ind){
          type_I_error <- type_I_error - exp(sapply(theta+delta, dbinom, x=boundary[ind[i]], size=n2, log=T)) * 
            exp(sapply(theta, dbinom, x=ind[i]-1, size=n1, log=T))
          boundary[ind[i]] <- boundary[ind[i]]-1
          if(boundary[ind[i]]<0){
            boundary_value[ind[i]] <- Inf
            if(min(boundary_value)==Inf){
              indicator <- 0
              break;
            }
          } else {
            boundary_value[ind[i]] <- func(ind[i]-1,boundary[ind[i]],n1,n2,delta,parameters)
          }
        }
        iteration <- iteration + 1
      }
    }
    if(indicator==0){
      print("Rejection region is the domain (reject the null all the time)")
      boundary <- rep(-1,n1+1)
      return(list(gamma, boundary, 0, iteration, boundary_initial, type_I_error_initial))
    } else{
      return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
    }
    
  } else{
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    boundary_value <- bayes_ts_value(boundary, n1, n2, delta, method, parameters)
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    while((max(type_I_error)>alpha)&indicator){
      ind_null <- which.max(type_I_error)
      while(type_I_error[ind_null]>alpha){
        gamma <- min(boundary_value)
        ind <- which(boundary_value==gamma) 
        n_ind <- length(ind)
        for(i in 1:n_ind){
          type_I_error <- type_I_error - exp(sapply(theta, dbinom, x=boundary[ind[i]], size=n1, log=T)) * 
            exp(sapply(theta+delta, dbinom, x=ind[i]-1, size=n2, log=T))
          boundary[ind[i]] <- boundary[ind[i]]+1
          if(boundary[ind[i]]>n1){
            boundary_value[ind[i]] <- Inf
            if(min(boundary_value)==Inf){
              indicator <- 0
              break
            }
          } else {
            boundary_value[ind[i]] <- func(boundary[ind[i]],ind[i]-1,n1,n2,delta,parameters)
          }
        }
        iteration <- iteration + 1
      }
    }
    if(indicator==0){
      print("Rejection region is the domain (reject the null all the time)")
      boundary <- rep(n1+1,n2+1)
      return(list(gamma, boundary, 0, iteration, boundary_initial, type_I_error_initial))
    } else{
      return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
    }
    
  }
  
}



##### Backward Procedure (Posterior probability) #####
backward_bayes <- function(n1, n2, delta, alpha, 
                           boundary_initial, type_I_error_initial, theta, 
                           method="Independent Beta", parameters=c(1,1,1,1)){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin with null: p1 <= p2 -delta
  #  -1 < delta < 1
  # alpha: intended significance level
  # boundary_initial (vector): boundary of the rejection region obtained based on the initial gamma
  # type_I_error_initial (vector): type I error at each grid based on the initial gamma
  # theta (vector): boundary of p1 with p2=p1+delta
  # method = "", and parameters of the corresponding method = parameters
  # Independent Beta (two independent beta prior with parameters a1, b1, a2, b2)
  # Logit Normal (logit-normal prior with parameters mu1, sigma1 (sd), mu2, sigma2(sd), sigma12 (cov))
  
  # select the test statistic
  if(method=="Independent Beta"){
    func <- bayes_indp_beta
  } else if(method=="Logit Normal"){
    func <- bayes_logit
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  # select which to update b_{y1}, or a_{y2}
  if( n2>n1 ){
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    candidate <- boundary_initial + 1
    boundary_value <- bayes_ts_value(candidate, n1, n2, delta, method, parameters)
    boundary_value[candidate>n2] <- -Inf
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    indicator <- 1 
    while(indicator==1){
      gamma_temp <- max(boundary_value)
      ind_temp <- which(boundary_value==gamma_temp)
      n_temp <- length(ind_temp)
      type_I_error_temp <- type_I_error
      
      # update the value required for next iteration
      for(i in 1:n_temp){
        type_I_error_temp <- type_I_error_temp + exp(sapply(theta+delta, dbinom, x=candidate[ind_temp[i]], size=n2, log=T)) * 
          exp(sapply(theta, dbinom, x=ind_temp[i]-1, size=n1, log=T)) 
        candidate[ind_temp[i]] <- candidate[ind_temp[i]]+1
        if(candidate[ind_temp[i]]>n2){
          boundary_value[ind_temp[i]] <- -Inf
        } else {
          boundary_value[ind_temp[i]] <- func(ind_temp[i]-1,candidate[ind_temp[i]],n1,n2,delta,parameters)
        }
      }
      
      if(max(type_I_error_temp)>alpha){
        indicator <- 0
        if(iteration==0){
          print("gamma does not change")
          return(list("No", boundary_initial, type_I_error_initial, 0, boundary_initial, type_I_error_initial))
        } else{
          gamma <- gamma_temp
          return(list(gamma, boundary, type_I_error, iteration, boundary_initial, type_I_error_initial))
        }
      } else{
        boundary[ind_temp] <- boundary[ind_temp] + 1
        type_I_error <- type_I_error_temp
      }
      iteration <- iteration + 1
    }
    
  } else{
    
    iteration <- 0
    
    # Initiate the value
    boundary <- boundary_initial
    candidate <- boundary_initial - 1
    boundary_value <- bayes_ts_value(candidate, n1, n2, delta, method, parameters)
    boundary_value[candidate<0] <- -Inf
    type_I_error <- type_I_error_initial
    
    # update the rejection region
    indicator <- 1
    while(indicator==1){
      gamma_temp <- max(boundary_value)
      ind_temp <- which(boundary_value==gamma_temp)
      n_temp <- length(ind_temp)
      type_I_error_temp <- type_I_error
      
      for(i in 1:n_temp){
        type_I_error_temp <- type_I_error_temp + exp(sapply(theta, dbinom, x=candidate[ind_temp[i]], size=n1, log=T)) * 
          exp(sapply(theta+delta, dbinom, x=ind_temp[i]-1, size=n2, log=T)) 
        candidate[ind_temp[i]] <- candidate[ind_temp[i]]-1
        if(candidate[ind_temp[i]]<0){
          boundary_value[ind_temp[i]] <- -Inf
        } else {
          boundary_value[ind_temp[i]] <- func(candidate[ind_temp[i]],ind_temp[i]-1,n1,n2,delta,parameters)
        }
      }
      
      if(max(type_I_error_temp)>alpha){
        indicator <- 0 
        if(iteration==0){
          print("gamma does not change")
          return(list("No", boundary_initial, type_I_error_initial, 0, boundary_initial, type_I_error_initial))
        } else{
          gamma <- gamma_temp
          return(list(gamma, boundary, type_I_error, iteration+1, boundary_initial, type_I_error_initial))
        }
      } else{
        boundary[ind_temp] <- boundary[ind_temp] - 1
        type_I_error <- type_I_error_temp
      }
      iteration <- iteration + 1
    }
    
  }
  
}




