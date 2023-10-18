######################################################################
###### Calculate the Type I error at the Boundary of Hypotheses ######
######################################################################
library(foreach)
library(doParallel)

typeIerror <- function(boundary, n1, n2, delta, theta){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # method: which statistic used
  # theta: a sequence of p1 at the boundary, with p2=p1+delta
  
  if( (delta>=0)&( (min(theta) < 0)||(max(theta) > 1-delta) )){
    print("The specifications of parameters are out of range!")
    return(NA)
  } 
  if( (delta<0)&( (min(theta) < -delta)||(max(theta) > 1) )){
    print("The specifications of parameters are out of range!")
    return(NA)
  }
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  if((length(boundary)==n1+1)&(n1!=n2)){
    # In this case, the boundary is for y2 given a y1
    return(foreach(i=theta, .combine=c) %dopar% {
      sum( sapply(seq(0,n1,1), dbinom, size=n1, prob=i)*sapply(boundary, pbinom, size=n2, prob=i+delta) )
    }) 
  } else{
    # In this case, the boundary is for y1 given a y2
    return(foreach(i=theta, .combine=c) %dopar% {
      sum( sapply(seq(0,n2,1), dbinom, size=n2, prob=i+delta)*(1-sapply(boundary-1, pbinom, size=n1, prob=i)) )
    })
  }
  
  stopImplicitCluster()
  
}



###########################################################################################
##### Calculate the Values of Test Statistics at the Boundary of the Rejection Region #####
###########################################################################################
### for frequentist method
# The out put is a (n2+1 x 1) vector or (n1+1 x 1) vector 
freq_ts_value <- function(boundary, n1, n2, delta, method="Unpooled Z"){
  
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
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  if(((length(boundary)==n1+1)&(n1!=n2))){
    return(foreach(i=0:n1, .combine=c) %dopar% {
      if(boundary[i+1]==-1){
        Inf
      } else{
        func(i, boundary[i+1], n1, n2, delta)
      }
    })
  } else{
    return(foreach(i=0:n2, .combine=c) %dopar% {
      if(boundary[i+1]==n1+1){
        Inf
      } else{
        func(boundary[i+1], i, n1, n2, delta)
      }
    })
  }
  
  stopImplicitCluster()
  
}

### for Bayesian method
# The out put is a (n2+1 x 1) vector or (n1+1 x 1) vector 
bayes_ts_value <- function(boundary, n1, n2, delta, method="Independent Beta", parameters=c(1,1,1,1)){
  
  if(method=="Independent Beta"){
    func <- bayes_indp_beta
  } else if(method=="Logit Normal"){
    func <- bayes_logit
  } else {
    print("This method is not available now!")
    return(NA)
  }
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  if(((length(boundary)==n1+1)&(n1!=n2))){
    return(foreach(i=0:n1, .combine=c) %dopar% {
      if(boundary[i+1]==-1){
        Inf
      } else{
        func(i, boundary[i+1], n1, n2, delta, parameters)
      }
    })
  } else{
    return(foreach(i=0:n2, .combine=c) %dopar% {
      if(boundary[i+1]==n1+1){
        Inf
      } else{
        func(boundary[i+1], i, n1, n2, delta, parameters)
      }
    })
  }
  
  stopImplicitCluster()
  
}



##################################################################
###### Calculate the Power of a Determinisitic Decision Rule #####
##################################################################
power_calculate <- function(boundary, n1, n2, delta, theta1, theta2){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # theta1, theta2: can be vectors for parameters
  # return will be a length(theta1)*length(theta2) matrix
  
  numCores <- detectCores()
  registerDoParallel(max(1, numCores-2))
  
  dim1 <- length(theta1)
  dim2 <- length(theta2)
  
  if((length(boundary)==n1+1)&(n1!=n2)){
    # In this case, the boundary is for y2 given a y1
    return(foreach(i=1:dim1, .combine=rbind) %dopar% {
      foreach(j=1:dim2, .combine=c)%dopar%{
        sum( sapply(seq(0,n1,1), dbinom, size=n1, prob=theta1[i])*sapply(boundary, pbinom, size=n2, prob=theta2[j]) )
      }
    })
  } else{
    # In this case, the boundary is for y1 given a y2
    return(foreach(i=1:dim1, .combine=rbind) %dopar% {
      foreach(j=1:dim2, .combine=c)%dopar%{
        sum( sapply(seq(0,n2,1), dbinom, size=n2, prob=theta2[j])*(1-sapply(boundary-1, pbinom, size=n1, prob=theta1[i])) )
      }
    })
  }
  
  stopImplicitCluster()
  
}



