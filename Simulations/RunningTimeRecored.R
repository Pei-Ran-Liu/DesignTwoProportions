#######################################
##### Running Time in Grid Search #####
#######################################

##### Superiority Trial #####
# Pooled Z
opc_zpool_sup <- c()
overalltime_zpool_sup <- c()
for(i in 1:5){
  overalltime_zpool_sup <- rbind(overalltime_zpool_sup, 
                                 system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method="Pooled Z",
                                                                    gamma=qnorm(0.975),seq(0.001,0.999,0.001),alpha=0.025,step=0.0001)))
  opc_zpool_sup <- rbind(opc_zpool_sup, res)
}
opc_zpoolalt_sup <- c()
overalltime_zpoolalt_sup <- c()
for(i in 1:5){
  overalltime_zpoolalt_sup <- rbind(overalltime_zpoolalt_sup, 
                                    system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method="Pooled Z",
                                                                       gamma=qnorm(0.975),seq(0.001,0.999,0.001),alpha=0.025,step=0.0001, zmethod=TRUE)))
  opc_zpoolalt_sup <- rbind(opc_zpoolalt_sup, res)
}

# Unpooled Z
opc_zunpool_sup <- c()
overalltime_zunpool_sup <- c()
for(i in 1:5){
  overalltime_zunpool_sup <- rbind(overalltime_zunpool_sup, 
                                   system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method="Unpooled Z",
                                                                      gamma=qnorm(0.975),seq(0.001,0.999,0.001),alpha=0.025,step=0.0001)))
  opc_zunpool_sup <- rbind(opc_zunpool_sup, res)
}
opc_zunpoolalt_sup <- c()
overalltime_zunpoolalt_sup <- c()
for(i in 1:5){
  overalltime_zunpoolalt_sup <- rbind(overalltime_zunpoolalt_sup, 
                                      system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method="Unpooled Z",
                                                                         gamma=qnorm(0.975),seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,zmethod=TRUE)))
  opc_zunpoolalt_sup <- rbind(opc_zunpoolalt_sup, res)
}

# Fisher's exact
opc_FisherCp_sup <- c()
overalltime_FisherCp_sup <- c()
for(i in 1:5){
  overalltime_FisherCp_sup <- rbind(overalltime_FisherCp_sup, 
                                    system.time(res <-girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method="Fisher's condition",
                                                                      gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001)))
  opc_FisherCp_sup <- rbind(opc_FisherCp_sup, res)
}

# Bayesian (Jeffreys)
opc_Jbayes_sup <- c()
overalltime_Jbayes_sup <- c()
for(i in 1:5){
  overalltime_Jbayes_sup <- rbind(overalltime_Jbayes_sup, 
                                  system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),
                                                                     gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,parameters=c(0.5,0.5,0.5,0.5))))
  opc_Jbayes_sup <- rbind(opc_Jbayes_sup, res)
}

# Bayes (Uniform) 
opc_Ubayes_sup <- c()
overalltime_Ubayes_sup <- c()
for(i in 1:5){
  overalltime_Ubayes_sup <- rbind(overalltime_Ubayes_sup, 
                                  system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),
                                                                     gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,parameters=c(1,1,1,1))))
  opc_Ubayes_sup <- rbind(opc_Ubayes_sup, res)
}

# Bayes (Logit N1) 
opc_Iln1_bayes_sup <- c()
overalltime_Iln1_bayes_sup <- c()
for(i in 1:5){
  overalltime_Iln1_bayes_sup <- rbind(overalltime_Iln1_bayes_sup, 
                                      system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                         gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,parameters=c(0,2,0,2,0))))
  opc_Iln1_bayes_sup <- rbind(opc_Iln1_bayes_sup, res)
}

# Bayes (Logit N2) 
opc_Iln2_bayes_sup <- c()
overalltime_Iln2_bayes_sup <- c()
for(i in 1:5){
  overalltime_Iln2_bayes_sup <- rbind(overalltime_Iln2_bayes_sup, 
                                      system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                         gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,parameters=c(0,3,0,3,0))))
  opc_Iln2_bayes_sup <- rbind(opc_Iln2_bayes_sup, res)
}

# Bayes (Logit NCorr) 
opc_Dln_bayes_sup <- c()
overalltime_Dln_bayes_sup <- c()
for(i in 1:5){
  overalltime_Dln_bayes_sup <- rbind(overalltime_Dln_bayes_sup, 
                                     system.time(res <- girds_algorithm(sup_design[i]/3*2,sup_design[i]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                        gamma=0.975,seq(0.001,0.999,0.001),alpha=0.025,step=0.0001,parameters=c(1,3,-1,3,0.5))))
  opc_Dln_bayes_sup <- rbind(opc_Dln_bayes_sup, res)
}

##### Noninferiority Trial #####
# Unpooled Z
opc_zunpool_nf <- c()
overalltime_zunpool_nf <- c()
for(i in 1:5){
  overalltime_zunpool_nf <- rbind(overalltime_zunpool_nf, 
                                  system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method="Unpooled Z",
                                                                     gamma=qnorm(0.975),seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001)))
  opc_zunpool_nf <- rbind(opc_zunpool_nf, res)
}
opc_zunpoolalt_nf <- c()
overalltime_zunpoolalt_nf <- c()
for(i in 1:5){
  overalltime_zunpoolalt_nf <- rbind(overalltime_zunpoolalt_nf, 
                                     system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method="Unpooled Z",
                                                                        gamma=qnorm(0.975),seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,zmethod=TRUE)))
  opc_zunpoolalt_nf <- rbind(opc_zunpoolalt_nf, res)
}

# Score
opc_score_nf <- c()
overalltime_score_nf <- c()
for(i in 1:5){
  overalltime_score_nf <- rbind(overalltime_score_nf, 
                                system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method="score",
                                                                   gamma=qnorm(0.975),seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001)))
  opc_score_nf <- rbind(opc_score_nf, res)
}
opc_scorealt_nf <- c()
overalltime_scorealt_nf <- c()
for(i in 1:5){
  overalltime_scorealt_nf <- rbind(overalltime_scorealt_nf, 
                                   system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method="score",
                                                                      gamma=qnorm(0.975),seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001, zmethod=TRUE)))
  opc_scorealt_nf <- rbind(opc_scorealt_nf, res)
}

# Bayesian (Jeffreys)
opc_Jbayes_nf <- c()
overalltime_Jbayes_nf <- c()
for(i in 1:5){
  overalltime_Jbayes_nf <- rbind(overalltime_Jbayes_nf, 
                                 system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),
                                                                    gamma=0.975,seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,parameters=c(0.5,0.5,0.5,0.5))))
  opc_Jbayes_nf <- rbind(opc_Jbayes_nf, res)
}

# Bayes (Uniform) 
opc_Ubayes_nf <- c()
overalltime_Ubayes_nf <- c()
for(i in 1:5){
  overalltime_Ubayes_nf <- rbind(overalltime_Ubayes_nf, 
                                 system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),
                                                                    gamma=0.975,seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,parameters=c(1,1,1,1))))
  opc_Ubayes_nf <- rbind(opc_Ubayes_nf, res)
}

# Bayes (Logit N1) 
opc_Iln1_bayes_nf <- c()
overalltime_Iln1_bayes_nf <- c()
for(i in 1:5){
  overalltime_Iln1_bayes_nf <- rbind(overalltime_Iln1_bayes_nf, 
                                     system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                        gamma=0.975,seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,parameters=c(0,2,0,2,0))))
  opc_Iln1_bayes_nf <- rbind(opc_Iln1_bayes_nf, res)
}

# Bayes (Logit N2) 
opc_Iln2_bayes_nf <- c()
overalltime_Iln2_bayes_nf <- c()
for(i in 1:5){
  overalltime_Iln2_bayes_nf <- rbind(overalltime_Iln2_bayes_nf, 
                                     system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                        gamma=0.975,seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,parameters=c(0,3,0,3,0))))
  opc_Iln2_bayes_nf <- rbind(opc_Iln2_bayes_nf, res)
}

# Bayes (Logit NCorr) 
opc_Dln_bayes_nf <- c()
overalltime_Dln_bayes_nf <- c()
for(i in 1:5){
  overalltime_Dln_bayes_nf <- rbind(overalltime_Dln_bayes_nf, 
                                    system.time(res <- girds_algorithm(nf_design[i]/4,nf_design[i]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),
                                                                       gamma=0.975,seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,step=0.0001,parameters=c(1,3,-1,3,0.5))))
  opc_Dln_bayes_nf <- rbind(opc_Dln_bayes_nf, res)
}


