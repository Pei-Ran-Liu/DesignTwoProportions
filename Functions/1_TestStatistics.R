#####################################################
##### Test Statistics in Two Sample Proportions #####
#####################################################

##### Pooled Z Statistic #####
z_pool <- function(y1, y2, n1, n2, delta){
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  
  p1 <- y1/n1
  p2 <- y2/n2
  pool <- (y1+y2)/(n1+n2)
  
  # for delta==0
  check_1 <- (delta==0) & ((y1==0&y2==0)||(y1==n1&y2==n2))
  # for delta!=0
  check_2 <- (delta!=0) & (y1==n1&y2==n2)
  
  if( check_1||check_2 ){
    return(0)
  } else {
    return( (p1-p2+delta)/sqrt(pool*(1-pool)*(1/n1+1/n2)) )
  }
  
}



##### Unpooled Z statistic #####
z_unpool <- function(y1, y2, n1, n2, delta){
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  
  p1 <- y1/n1
  p2 <- y2/n2
  
  # for delta==0
  check_1 <- (delta==0) & (y1==0&y2==0||y1==n1&y2==n2)
  # for another delta==0
  check_2 <- (delta==0) & (y1==n1&y2==0||y1==0&y2==n2)
  # for delta!=0
  check_3 <- (delta!=0) & (y1==n1&y2==n2)
  
  if( check_1||check_3 ){
    return(0)
  } else if(check_2){
    if(y1==n1&y2==0){
      return(Inf)
    } else{
      return(-Inf)
    }
  } else {
    return( (p1-p2+delta)/sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2) )
  }
  
}



##### Score Statistic #####
### Boundary case ###
y1_0 <- function(y1, y2, n1, n2, delta){
  delta0 <- ( (n1+n2)-sqrt((n1+n2)^2-4*n1*y2) ) / (2*n1)
  if(delta>=delta0){
    return(0)
  } else{
    a <- n1+n2
    b <- 2*delta*n1-n1-y2-n2+delta*n2
    c <- delta^2*n1-delta*(n1+n2)+y2
    return( (-b-sqrt(b^2-4*a*c))/(2*a) )
  }
}
left_uppercorner <- function(y1, y2, n1, n2, delta){
  if(y1==0){
    return(y1_0(y1, y2, n1, n2, delta))
  } else if(y2==n2){
    return(1-y1_0(n2-y2, n1-y1, n2, n1, delta)-delta) # y1_0 function get 1-tilde(p_2)
  } else{
    print("You are in a wrong function!")
  }
}
y1_n1 <- function(y1, y2, n1, n2, delta){
  delta0 <- -( (n1+n2)-sqrt((n1+n2)^2+4*n1*(y2-n2)) ) / (2*n1)
  if(delta<=delta0){
    return(1)
  } else{
    a <- -(n1+n2)
    b <- n1-2*delta*n1+y2-delta*n2
    c <- -delta^2*n1+delta*n1
    return( (-b-sqrt(b^2-4*a*c))/(2*a) )
  }
}
right_lowercorner <- function(y1, y2, n1, n2, delta){
  if(y1==n1){
    return(y1_n1(y1, y2, n1, n2, delta))
  } else if(y2==0){
    return(1-y1_n1(n2-y2, n1-y1, n2, n1, delta)-delta)
  } else{
    print("You are in a wrong function!")
  }
} 

### Overall ###
score <- function(y1, y2, n1, n2, delta){
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  
  p1 <- y1/n1
  p2 <- y2/n2
  
  # whether it is at the boundary
  boundary_check <- (y1==0)||(y1==n1)||(y2==0)||(y2==n2)
  
  if(delta==0){
    # when margin == 0, it turns the pooled z test statistic
    return(z_pool(y1,y2,n1,n2,delta))
  } else if(!boundary_check){
    # when none of observations at the boundary, it retures the closed form solution
    theta <- n2/n1
    a <- 1+theta
    b <- -(1+theta+p1+theta*p2-delta*(theta+2))
    c <- delta^2-delta*(2*p1+theta+1)+p1+theta*p2
    d <- p1*delta*(1-delta)
    v <- b^3/((3*a)^3)-b*c/(6*a^2)+d/(2*a)
    u <- sign(v)*sqrt(b^2/((3*a)^2)-c/(3*a))
    if(u==0){
      p1_rmle <- -b/(3*a)
    } else{
      w <- (pi+acos(v/(u^3)))/3
      p1_rmle <- 2*u*cos(w)-b/(3*a)
    }
    p2_rmle <- p1_rmle + delta
    return( (p1-p2+delta)/sqrt(p1_rmle*(1-p1_rmle)/n1+p2_rmle*(1-p2_rmle)/n2) )
  } else{
    # boundary cases
    if((y1==0)&(y2==n2)){
      return(-Inf)
    } else if((y1==n1)&(y2==0)){
      return(Inf)
    } else if((y1==0)&(y2==0)){
      p1_rmle <- max(0, -delta)
      p2_rmle <- p1_rmle + delta
      return( (p1-p2+delta)/sqrt(p1_rmle*(1-p1_rmle)/n1+p2_rmle*(1-p2_rmle)/n2) )
    } else if((y1==n1)&(y2==n2)){
      p1_rmle <- min(1,1-delta)
      p2_rmle <- p1_rmle + delta
      return( (p1-p2+delta)/sqrt(p1_rmle*(1-p1_rmle)/n1+p2_rmle*(1-p2_rmle)/n2) )
    } else if((y1==0)||(y2==n2)){
      p1_rmle <- left_uppercorner(y1, y2, n1, n2, delta)
      p2_rmle <- p1_rmle + delta
      return( (p1-p2+delta)/sqrt(p1_rmle*(1-p1_rmle)/n1+p2_rmle*(1-p2_rmle)/n2) )
    } else if ((y1==n1)||(y2==0)){
      p1_rmle <- right_lowercorner(y1, y2, n1, n2, delta)
      p2_rmle <- p1_rmle + delta
      return( (p1-p2+delta)/sqrt(p1_rmle*(1-p1_rmle)/n1+p2_rmle*(1-p2_rmle)/n2) )
    } 
  } 
}



##### Fisher Conditional p value #####
fisher_condp <- function(y1, y2, n1, n2, delta){
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  
  if(delta==0){
    return( exp(phyper(y1-1, n1, n2, y1+y2, log.p=T)) )
  } else{
    print("Fisher's exact test does not applicable for delta \ne 0.")
    return(NA)
  }
}



##### Bayesian Method #####
library(distrEx)
# Independent Beta Prior
bayes_indp_beta <- function(y1, y2, n1, n2, delta, parameters){
  
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  # parameters
  # a1, b1: beta prior parameter for the population 1
  # a2, b2: beta prior parameter for the population 2
  
  a1 <- parameters[1]
  b1 <- parameters[2]
  a2 <- parameters[3]
  b2 <- parameters[4]
  
  integrand <- function(p2){
    exp(dbeta(p2, y2+a2, n2-y2+b2, log = T) + pbeta(p2-delta, y1+a1, n1-y1+b1, lower.tail = F, log.p = T))
  }
  
  if(delta>0){
    part1 <- GLIntegrate(integrand,delta,1)
    part2 <- pbeta(delta, y2+a2, n2-y2+b2)
    return(part1+part2)
  } else{
    return(GLIntegrate(integrand,0,1+delta))
  }
  
}

# Logit Normal Prior (turns to independent prior given sigma12=0)
bayes_logit <- function(y1, y2, n1, n2, delta, parameters){
  
  # y1: response for the population 1
  # y2: response for the population 2
  # n1: total number of the population 1
  # n2: total number of the population 2
  # delta: test margin for the null: p1 <= p2 - delta
  # -1 < delta < 1
  # parameters
  # mu1, sigma1: normal prior parameter for the population 1 on a logit scale (sd)
  # mu2, sigma2: normal prior parameter for the population 1 on a logit scale (sd)
  # sigma12: covariance between two parameters
  
  mu1 <- parameters[1]
  sigma1 <- parameters[2]
  mu2 <- parameters[3]
  sigma2 <- parameters[4]
  sigma12 <- parameters[5]
  
  ### numernator
  part_p1 <- function(p2){
    sapply(p2, function(z) {GLIntegrate(function(p1) {
      exp(-(log(p1/(1-p1))-(mu1+sigma12*sigma2^(-2)*(z-mu2)))^2/2/(sigma1^2-sigma12^2/sigma2^2)+
            (y1-1)*log(p1)+(n1-y1-1)*log(1-p1))
    }, z-delta, 1)*
        exp(-(log(z/(1-z))-mu2)^2/2/sigma2^2+(y2-1)*log(z)+(n2-y2-1)*log(1-z))
    })
  }
  part_p2 <- function(p2){
    sapply(p2, function(z) {GLIntegrate(function(p1) {
      exp(-(log(p1/(1-p1))-(mu1+sigma12*sigma2^(-2)*(z-mu2)))^2/2/(sigma1^2-sigma12^2/sigma2^2)+
            (y1-1)*log(p1)+(n1-y1-1)*log(1-p1))
    }, 0, 1)*
        exp(-(log(z/(1-z))-mu2)^2/2/sigma2^2+(y2-1)*log(z)+(n2-y2-1)*log(1-z))
    })
  }
  
  if(delta>0){
    part1 <- GLIntegrate(part_p1,delta,1)
    part2 <- GLIntegrate(part_p2,0,delta)
    num <- log(part1+part2)
  } else{
    num <- log(GLIntegrate(part_p1,0,1+delta))
  }
  
  ### denominator
  den <- log(GLIntegrate(part_p2,0,1))
  
  ### final output
  if(exp(num)==0){
    return(0)
  } else if(exp(num-den)>1){
    return(1)
  } else{
    return(exp(num-den))
  }
  
}


