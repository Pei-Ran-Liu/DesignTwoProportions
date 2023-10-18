############################################################
##### Combined Single Function for the Proposed Method #####
############################################################

calibration_algorithm <- function(n1, n2, delta, method="Unpooled Z", 
                                  gamma_initial=0.95, theta, alpha=0.05, 
                                  parameters=list(c(1,1,1,1),c(1,1,1,1))){
  
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin with null: p1 <= p2 - delta
  # method: which statistic used in a test or a design 
  #         for freqeuntist method, indicate the method directly
  #         for Bayesian posterior probability method, the second component indicate the method
  # gamma_initial: initial cut value of a decision rule
  # theta: a sequence of p1 at the boundary, with p2=p1+delta
  # Note: this theta is the range that you want to control the size rather than to calculate the exact p value
  # Also, the range of theta is important
  # delta >=0: theta \in (0,1-delta)
  # delta <=0: theta \in (-delta,1)
  # alpha: intended significance level
  
  if(method[1]=="Bayes Posterior Probability"){
    
    # Initialize the rejection region and calculate the size 
    rr_initial <- bayes_rr(n1, n2, delta, gamma_initial, method[2], parameters)
    type_I_error_initial <- typeIerror(rr_initial, n1, n2, delta, theta)
    max_initial <- max(type_I_error_initial)
    
    # Select a proper procedure and tuning
    if(max_initial==alpha){
      return(list(gamma_initial, rr_initial, type_I_error_initial, 0, rr_initial, type_I_error_initial))
    } else if(max_initial>alpha){
      return(forward_bayes(n1, n2, delta, alpha, rr_initial, type_I_error_initial, theta, method=method[2], parameters))
    } else{
      res <- backward_bayes(n1, n2, delta, alpha, rr_initial, type_I_error_initial, theta, method=method[2], parameters)
      if(res[[1]]!="No"){
        return(res)
      } else{
        res[[1]] <- gamma_initial
        return(res)
      }
    }
    
  } else{
    
    # Initialize the rejection region and calculate the size 
    rr_initial <- freq_rr(n1, n2, delta, gamma_initial, method)
    type_I_error_initial <- typeIerror(rr_initial, n1, n2, delta, theta)
    max_initial <- max(type_I_error_initial)
    
    # Select a proper procedure and tuning
    if(max_initial==alpha){
      return(list(gamma_initial, rr_initial, type_I_error_initial, 0, rr_initial,type_I_error_initial))
    } else if(max_initial>alpha){
      return(forward_freq(n1, n2, delta, alpha, rr_initial, type_I_error_initial, theta, method))
    } else{
      res <- backward_freq(n1, n2, delta, alpha, rr_initial, type_I_error_initial, theta, method)
      if(res[[1]]!="No"){
        return(res)
      } else{
        res[[1]] <- gamma_initial
        return(res)
      }
    }
    
  }
  
}


