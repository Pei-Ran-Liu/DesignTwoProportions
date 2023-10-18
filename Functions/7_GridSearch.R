################################################################
##### Grid Search Algorithm for Comparison in Running Time #####
################################################################

### evaluate the type I error for the current gamma
gs_evaluate <- function(n1, n2, delta, method="Unpooled Z", 
                        gamma, theta, alpha=0.05, parameters=list(c(1,1,1,1),c(1,1,1,1))){
  
  if(method[1]=="Bayes Posterior Probability"){
    
    rr <- bayes_rr(n1, n2, delta, gamma, method[2], parameters)
    type_I_error <- typeIerror(rr, n1, n2, delta, theta)
    return(max(type_I_error))
    
  } else{
    
    rr <- freq_rr(n1, n2, delta, gamma, method)
    type_I_error <- typeIerror(rr, n1, n2, delta, theta)
    return(max(type_I_error))
    
  }
  
}

### generate the final gamma based on grid search algorithm
girds_algorithm <- function(n1, n2, delta, method="Unpooled Z", gamma=0.95, theta, alpha, 
                            parameters=list(c(1,1,1,1),c(1,1,1,1)), step=0.0001,
                            zmethod=FALSE){
  initial_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
  
  if((method[1]%in%c("Unpooled Z", "Pooled Z", "score"))&(zmethod==TRUE)){
    
    if(initial_tIe>0.025){
      curret_tIe <- initial_tIe
      while(curret_tIe>0.025){
        gamma <- pnorm(gamma)
        gamma <- gamma + step
        gamma <- qnorm(gamma)
        curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
      }
    } else{
      curret_tIe <- initial_tIe
      while(curret_tIe<=0.025){
        gamma <- pnorm(gamma)
        gamma <- gamma - step
        gamma <- qnorm(gamma)
        curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
      }
      gamma <- gamma + step
      curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
    }
    return(c(gamma, curret_tIe))
    
  } else{
    
    if(initial_tIe>0.025){
      curret_tIe <- initial_tIe
      while(curret_tIe>0.025){
        gamma <- gamma + step
        curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
      }
    } else{
      curret_tIe <- initial_tIe
      while(curret_tIe<=0.025){
        gamma <- gamma - step
        curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
      }
      gamma <- gamma + step
      curret_tIe <- gs_evaluate(n1, n2, delta, method, gamma, theta, alpha, parameters)
    }
    return(c(gamma, curret_tIe))
    
  }
  
}


