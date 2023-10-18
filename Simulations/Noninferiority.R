##############################################
##### Design of the Noninferiority Trial #####
##############################################

# Intended power: 0.9 at 0.06, 0.06; Allocation ratio 1:3
nf_design <- seq(1320, 1640, 80)

##### Unpooled Z #####
zunpool_nf_time_overall <- c()
zunpool_nf_time_overall <- rbind(zunpool_nf_time_overall,
                                 system.time(zunpool_design1 <- 
                                               calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Unpooled Z",gamma_initial=qnorm(0.95),
                                                                     theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
zunpool_nf_time_overall <- rbind(zunpool_nf_time_overall,
                                 system.time(zunpool_design2 <- 
                                               calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Unpooled Z",gamma_initial=qnorm(0.95),
                                                                     theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
zunpool_nf_time_overall <- rbind(zunpool_nf_time_overall,
                                 system.time(zunpool_design3 <- 
                                               calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Unpooled Z",gamma_initial=qnorm(0.95),
                                                                     theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
zunpool_nf_time_overall <- rbind(zunpool_nf_time_overall,
                                 system.time(zunpool_design4 <- 
                                               calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Unpooled Z",gamma_initial=qnorm(0.95),
                                                                     theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
zunpool_nf_time_overall <- rbind(zunpool_nf_time_overall,
                                 system.time(zunpool_design5 <- 
                                               calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Unpooled Z",gamma_initial=qnorm(0.95),
                                                                     theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
# randomized part
zunpool_design1_rr_part2 <- secondp_rd_freq(zunpool_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Unpooled Z")
zunpool_design1_tIe_part2 <- powerf_2nd(zunpool_design1[[2]], zunpool_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
zunpool_design1_xi <- binarys_xi(zunpool_design1[[3]], zunpool_design1_tIe_part2, 0.05, precision = 0.00001)

zunpool_design2_rr_part2 <- secondp_rd_freq(zunpool_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Unpooled Z")
zunpool_design2_tIe_part2 <- powerf_2nd(zunpool_design2[[2]], zunpool_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
zunpool_design2_xi <- binarys_xi(zunpool_design2[[3]], zunpool_design2_tIe_part2, 0.05, precision = 0.00001)

zunpool_design3_rr_part2 <- secondp_rd_freq(zunpool_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Unpooled Z")
zunpool_design3_tIe_part2 <- powerf_2nd(zunpool_design3[[2]], zunpool_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
zunpool_design3_xi <- binarys_xi(zunpool_design3[[3]], zunpool_design3_tIe_part2, 0.05, precision = 0.00001)

zunpool_design4_rr_part2 <- secondp_rd_freq(zunpool_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Unpooled Z")
zunpool_design4_tIe_part2 <- powerf_2nd(zunpool_design4[[2]], zunpool_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
zunpool_design4_xi <- binarys_xi(zunpool_design4[[3]], zunpool_design4_tIe_part2, 0.05, precision = 0.00001)

zunpool_design5_rr_part2 <- secondp_rd_freq(zunpool_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Unpooled Z")
zunpool_design5_tIe_part2 <- powerf_2nd(zunpool_design5[[2]], zunpool_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
zunpool_design5_xi <- binarys_xi(zunpool_design5[[3]], zunpool_design5_tIe_part2, 0.05, precision = 0.00001)



##### Score #####
score_nf_time_overall <- c()
score_nf_time_overall <- rbind(score_nf_time_overall,
                               system.time(score_design1 <- 
                                             calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="score",gamma_initial=qnorm(0.95),
                                                                   theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
score_nf_time_overall <- rbind(score_nf_time_overall,
                               system.time(score_design2 <- 
                                             calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="score",gamma_initial=qnorm(0.95),
                                                                   theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
score_nf_time_overall <- rbind(score_nf_time_overall,
                               system.time(score_design3 <- 
                                             calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="score",gamma_initial=qnorm(0.95),
                                                                   theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
score_nf_time_overall <- rbind(score_nf_time_overall,
                               system.time(score_design4 <- 
                                             calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="score",gamma_initial=qnorm(0.95),
                                                                   theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
score_nf_time_overall <- rbind(score_nf_time_overall,
                               system.time(score_design5 <- 
                                             calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="score",gamma_initial=qnorm(0.95),
                                                                   theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05)))
# randomized part
score_design1_rr_part2 <- secondp_rd_freq(score_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="score")
score_design1_tIe_part2 <- powerf_2nd(score_design1[[2]], score_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
score_design1_xi <- binarys_xi(score_design1[[3]], score_design1_tIe_part2, 0.05, precision = 0.00001)

score_design2_rr_part2 <- secondp_rd_freq(score_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="score")
score_design2_tIe_part2 <- powerf_2nd(score_design2[[2]], score_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
score_design2_xi <- binarys_xi(score_design2[[3]], score_design2_tIe_part2, 0.05, precision = 0.00001)

score_design3_rr_part2 <- secondp_rd_freq(score_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="score")
score_design3_tIe_part2 <- powerf_2nd(score_design3[[2]], score_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
score_design3_xi <- binarys_xi(score_design3[[3]], score_design3_tIe_part2, 0.05, precision = 0.00001)

score_design4_rr_part2 <- secondp_rd_freq(score_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="score")
score_design4_tIe_part2 <- powerf_2nd(score_design4[[2]], score_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
score_design4_xi <- binarys_xi(score_design4[[3]], score_design4_tIe_part2, 0.05, precision = 0.00001)

score_design5_rr_part2 <- secondp_rd_freq(score_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="score")
score_design5_tIe_part2 <- powerf_2nd(score_design5[[2]], score_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
score_design5_xi <- binarys_xi(score_design5[[3]], score_design5_tIe_part2, 0.05, precision = 0.00001)



##### Bayesian (Jeffreys) #####
Jbayes_nf_time_overall <- c()
Jbayes_nf_time_overall <- rbind(Jbayes_nf_time_overall,
                                system.time(Jbayes_design1 <- 
                                              calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_nf_time_overall <- rbind(Jbayes_nf_time_overall,
                                system.time(Jbayes_design2 <- 
                                              calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_nf_time_overall <- rbind(Jbayes_nf_time_overall,
                                system.time(Jbayes_design3 <- 
                                              calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_nf_time_overall <- rbind(Jbayes_nf_time_overall,
                                system.time(Jbayes_design4 <- 
                                              calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_nf_time_overall <- rbind(Jbayes_nf_time_overall,
                                system.time(Jbayes_design5 <- 
                                              calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0.5,0.5,0.5,0.5))))
# randomized part
Jbayes_design1_rr_part2 <- secondp_rd_bayes(Jbayes_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design1_tIe_part2 <- powerf_2nd(Jbayes_design1[[2]], Jbayes_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Jbayes_design1_xi <- binarys_xi(Jbayes_design1[[3]], Jbayes_design1_tIe_part2, 0.05, precision = 0.00001)

Jbayes_design2_rr_part2 <- secondp_rd_bayes(Jbayes_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design2_tIe_part2 <- powerf_2nd(Jbayes_design2[[2]], Jbayes_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Jbayes_design2_xi <- binarys_xi(Jbayes_design2[[3]], Jbayes_design2_tIe_part2, 0.05, precision = 0.00001)

Jbayes_design3_rr_part2 <- secondp_rd_bayes(Jbayes_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design3_tIe_part2 <- powerf_2nd(Jbayes_design3[[2]], Jbayes_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Jbayes_design3_xi <- binarys_xi(Jbayes_design3[[3]], Jbayes_design3_tIe_part2, 0.05, precision = 0.00001)

Jbayes_design4_rr_part2 <- secondp_rd_bayes(Jbayes_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design4_tIe_part2 <- powerf_2nd(Jbayes_design4[[2]], Jbayes_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Jbayes_design4_xi <- binarys_xi(Jbayes_design4[[3]], Jbayes_design4_tIe_part2, 0.05, precision = 0.00001)

Jbayes_design5_rr_part2 <- secondp_rd_bayes(Jbayes_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design5_tIe_part2 <- powerf_2nd(Jbayes_design5[[2]], Jbayes_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Jbayes_design5_xi <- binarys_xi(Jbayes_design5[[3]], Jbayes_design5_tIe_part2, 0.05, precision = 0.00001)



##### Bayes (Uniform) #####
Ubayes_nf_time_overall <- c()
Ubayes_nf_time_overall <- rbind(Ubayes_nf_time_overall,
                                system.time(Ubayes_design1 <- 
                                              calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,1,1,1))))
Ubayes_nf_time_overall <- rbind(Ubayes_nf_time_overall,
                                system.time(Ubayes_design2 <- 
                                              calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,1,1,1))))
Ubayes_nf_time_overall <- rbind(Ubayes_nf_time_overall,
                                system.time(Ubayes_design3 <- 
                                              calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,1,1,1))))
Ubayes_nf_time_overall <- rbind(Ubayes_nf_time_overall,
                                system.time(Ubayes_design4 <- 
                                              calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,1,1,1))))
Ubayes_nf_time_overall <- rbind(Ubayes_nf_time_overall,
                                system.time(Ubayes_design5 <- 
                                              calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.95,
                                                                    theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,1,1,1))))
# randomized part
Ubayes_design1_rr_part2 <- secondp_rd_bayes(Ubayes_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design1_tIe_part2 <- powerf_2nd(Ubayes_design1[[2]], Ubayes_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Ubayes_design1_xi <- binarys_xi(Ubayes_design1[[3]], Ubayes_design1_tIe_part2, 0.05, precision = 0.00001)

Ubayes_design2_rr_part2 <- secondp_rd_bayes(Ubayes_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design2_tIe_part2 <- powerf_2nd(Ubayes_design2[[2]], Ubayes_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Ubayes_design2_xi <- binarys_xi(Ubayes_design2[[3]], Ubayes_design2_tIe_part2, 0.05, precision = 0.00001)

Ubayes_design3_rr_part2 <- secondp_rd_bayes(Ubayes_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design3_tIe_part2 <- powerf_2nd(Ubayes_design3[[2]], Ubayes_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Ubayes_design3_xi <- binarys_xi(Ubayes_design3[[3]], Ubayes_design3_tIe_part2, 0.05, precision = 0.00001)

Ubayes_design4_rr_part2 <- secondp_rd_bayes(Ubayes_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design4_tIe_part2 <- powerf_2nd(Ubayes_design4[[2]], Ubayes_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Ubayes_design4_xi <- binarys_xi(Ubayes_design4[[3]], Ubayes_design4_tIe_part2, 0.05, precision = 0.00001)

Ubayes_design5_rr_part2 <- secondp_rd_bayes(Ubayes_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design5_tIe_part2 <- powerf_2nd(Ubayes_design5[[2]], Ubayes_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Ubayes_design5_xi <- binarys_xi(Ubayes_design5[[3]], Ubayes_design5_tIe_part2, 0.05, precision = 0.00001)



##### Bayes (Logit N1) #####
Iln1_bayes_nf_time_overall <- c()
Iln1_bayes_nf_time_overall <- rbind(Iln1_bayes_nf_time_overall,
                                    system.time(Iln1_bayes_design1 <- 
                                                  calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,2,0,2,0))))
Iln1_bayes_nf_time_overall <- rbind(Iln1_bayes_nf_time_overall,
                                    system.time(Iln1_bayes_design2 <- 
                                                  calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,2,0,2,0))))
Iln1_bayes_nf_time_overall <- rbind(Iln1_bayes_nf_time_overall,
                                    system.time(Iln1_bayes_design3 <- 
                                                  calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,2,0,2,0))))
Iln1_bayes_nf_time_overall <- rbind(Iln1_bayes_nf_time_overall,
                                    system.time(Iln1_bayes_design4 <- 
                                                  calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,2,0,2,0))))
Iln1_bayes_nf_time_overall <- rbind(Iln1_bayes_nf_time_overall,
                                    system.time(Iln1_bayes_design5 <- 
                                                  calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,2,0,2,0))))
# randomized part
Iln1_bayes_design1_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design1_tIe_part2 <- powerf_2nd(Iln1_bayes_design1[[2]], Iln1_bayes_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln1_bayes_design1_xi <- binarys_xi(Iln1_bayes_design1[[3]], Iln1_bayes_design1_tIe_part2, 0.05, precision = 0.00001)

Iln1_bayes_design2_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design2_tIe_part2 <- powerf_2nd(Iln1_bayes_design2[[2]], Iln1_bayes_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln1_bayes_design2_xi <- binarys_xi(Iln1_bayes_design2[[3]], Iln1_bayes_design2_tIe_part2, 0.05, precision = 0.00001)

Iln1_bayes_design3_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design3_tIe_part2 <- powerf_2nd(Iln1_bayes_design3[[2]], Iln1_bayes_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln1_bayes_design3_xi <- binarys_xi(Iln1_bayes_design3[[3]], Iln1_bayes_design3_tIe_part2, 0.05, precision = 0.00001)

Iln1_bayes_design4_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design4_tIe_part2 <- powerf_2nd(Iln1_bayes_design4[[2]], Iln1_bayes_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln1_bayes_design4_xi <- binarys_xi(Iln1_bayes_design4[[3]], Iln1_bayes_design4_tIe_part2, 0.05, precision = 0.00001)

Iln1_bayes_design5_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design5_tIe_part2 <- powerf_2nd(Iln1_bayes_design5[[2]], Iln1_bayes_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln1_bayes_design5_xi <- binarys_xi(Iln1_bayes_design5[[3]], Iln1_bayes_design5_tIe_part2, 0.05, precision = 0.00001)



##### Bayes (Logit N2) #####
Iln2_bayes_nf_time_overall <- c()
Iln2_bayes_nf_time_overall <- rbind(Iln2_bayes_nf_time_overall,
                                    system.time(Iln2_bayes_design1 <- 
                                                  calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,3,0,3,0))))
Iln2_bayes_nf_time_overall <- rbind(Iln2_bayes_nf_time_overall,
                                    system.time(Iln2_bayes_design2 <- 
                                                  calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,3,0,3,0))))
end_time <- Sys.time()
Iln2_bayes_nf_time <- c(Iln2_bayes_nf_time, difftime(end_time, start_time, units="secs"))

start_time <- Sys.time()
Iln2_bayes_nf_time_overall <- rbind(Iln2_bayes_nf_time_overall,
                                    system.time(Iln2_bayes_design3 <- 
                                                  calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,3,0,3,0))))
Iln2_bayes_nf_time_overall <- rbind(Iln2_bayes_nf_time_overall,
                                    system.time(Iln2_bayes_design4 <- 
                                                  calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,3,0,3,0))))
Iln2_bayes_nf_time_overall <- rbind(Iln2_bayes_nf_time_overall,
                                    system.time(Iln2_bayes_design5 <- 
                                                  calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                        theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(0,3,0,3,0))))
# randomized part
Iln2_bayes_design1_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design1_tIe_part2 <- powerf_2nd(Iln2_bayes_design1[[2]], Iln2_bayes_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln2_bayes_design1_xi <- binarys_xi(Iln2_bayes_design1[[3]], Iln2_bayes_design1_tIe_part2, 0.05, precision = 0.00001)

Iln2_bayes_design2_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design2_tIe_part2 <- powerf_2nd(Iln2_bayes_design2[[2]], Iln2_bayes_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln2_bayes_design2_xi <- binarys_xi(Iln2_bayes_design2[[3]], Iln2_bayes_design2_tIe_part2, 0.05, precision = 0.00001)

Iln2_bayes_design3_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design3_tIe_part2 <- powerf_2nd(Iln2_bayes_design3[[2]], Iln2_bayes_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln2_bayes_design3_xi <- binarys_xi(Iln2_bayes_design3[[3]], Iln2_bayes_design3_tIe_part2, 0.05, precision = 0.00001)

Iln2_bayes_design4_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design4_tIe_part2 <- powerf_2nd(Iln2_bayes_design4[[2]], Iln2_bayes_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln2_bayes_design4_xi <- binarys_xi(Iln2_bayes_design4[[3]], Iln2_bayes_design4_tIe_part2, 0.05, precision = 0.00001)

Iln2_bayes_design5_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design5_tIe_part2 <- powerf_2nd(Iln2_bayes_design5[[2]], Iln2_bayes_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Iln2_bayes_design5_xi <- binarys_xi(Iln2_bayes_design5[[3]], Iln2_bayes_design5_tIe_part2, 0.05, precision = 0.00001)



##### Bayes (Logit NCorr) #####
Dln_bayes_nf_time_overall <- c()
Dln_bayes_nf_time_overall <- rbind(Dln_bayes_nf_time_overall,
                                   system.time(Dln_bayes_design1 <- 
                                                 calibration_algorithm(nf_design[1]/4,nf_design[1]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                       theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_nf_time_overall <- rbind(Dln_bayes_nf_time_overall,
                                   system.time(Dln_bayes_design2 <- 
                                                 calibration_algorithm(nf_design[2]/4,nf_design[2]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                       theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_nf_time_overall <- rbind(Dln_bayes_nf_time_overall,
                                   system.time(Dln_bayes_design3 <- 
                                                 calibration_algorithm(nf_design[3]/4,nf_design[3]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                       theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_nf_time_overall <- rbind(Dln_bayes_nf_time_overall,
                                   system.time(Dln_bayes_design4 <- 
                                                 calibration_algorithm(nf_design[4]/4,nf_design[4]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                       theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_nf_time_overall <- rbind(Dln_bayes_nf_time_overall,
                                   system.time(Dln_bayes_design5 <- 
                                                 calibration_algorithm(nf_design[5]/4,nf_design[5]/4*3,4.1/100,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.95,
                                                                       theta=seq(0.001,1-4.1/100-0.001,0.001),alpha=0.05,parameters=c(1,3,-1,3,0.5))))
# randomized part
Dln_bayes_design1_rr_part2 <- secondp_rd_bayes(Dln_bayes_design1[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design1_tIe_part2 <- powerf_2nd(Dln_bayes_design1[[2]], Dln_bayes_design1_rr_part2, nf_design[1]/4,nf_design[1]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Dln_bayes_design1_xi <- binarys_xi(Dln_bayes_design1[[3]], Dln_bayes_design1_tIe_part2, 0.05, precision = 0.00001)

Dln_bayes_design2_rr_part2 <- secondp_rd_bayes(Dln_bayes_design2[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design2_tIe_part2 <- powerf_2nd(Dln_bayes_design2[[2]], Dln_bayes_design2_rr_part2, nf_design[2]/4,nf_design[2]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Dln_bayes_design2_xi <- binarys_xi(Dln_bayes_design2[[3]], Dln_bayes_design2_tIe_part2, 0.05, precision = 0.00001)

Dln_bayes_design3_rr_part2 <- secondp_rd_bayes(Dln_bayes_design3[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design3_tIe_part2 <- powerf_2nd(Dln_bayes_design3[[2]], Dln_bayes_design3_rr_part2, nf_design[3]/4,nf_design[3]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Dln_bayes_design3_xi <- binarys_xi(Dln_bayes_design3[[3]], Dln_bayes_design3_tIe_part2, 0.05, precision = 0.00001)

Dln_bayes_design4_rr_part2 <- secondp_rd_bayes(Dln_bayes_design4[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design4_tIe_part2 <- powerf_2nd(Dln_bayes_design4[[2]], Dln_bayes_design4_rr_part2, nf_design[4]/4,nf_design[4]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Dln_bayes_design4_xi <- binarys_xi(Dln_bayes_design4[[3]], Dln_bayes_design4_tIe_part2, 0.05, precision = 0.00001)

Dln_bayes_design5_rr_part2 <- secondp_rd_bayes(Dln_bayes_design5[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design5_tIe_part2 <- powerf_2nd(Dln_bayes_design5[[2]], Dln_bayes_design5_rr_part2, nf_design[5]/4,nf_design[5]/4*3,4.1/100, seq(0.001,1-4.1/100-0.001,0.001))
Dln_bayes_design5_xi <- binarys_xi(Dln_bayes_design5[[3]], Dln_bayes_design5_tIe_part2, 0.05, precision = 0.00001)



##### Operating Characteristics #####
### Type I Error ###
# unpooled z
sizenf_zunpool_before <- c(max(zunpool_design1_nf[[6]]),max(zunpool_design2_nf[[6]]),
                           max(zunpool_design3_nf[[6]]),max(zunpool_design4_nf[[6]]),
                           max(zunpool_design5_nf[[6]]))
sizenf_zunpool_after <- c(max(zunpool_design1_nf[[3]]),max(zunpool_design2_nf[[3]]),
                          max(zunpool_design3_nf[[3]]),max(zunpool_design4_nf[[3]]),
                          max(zunpool_design5_nf[[3]]))
sizenf_zunpool_random <- c(zunpool_design1_nf_xi[2], zunpool_design2_nf_xi[2],
                           zunpool_design3_nf_xi[2], zunpool_design4_nf_xi[2],
                           zunpool_design5_nf_xi[2])
# score
sizenf_score_before <- c(max(score_design1_nf[[6]]),max(score_design2_nf[[6]]),
                         max(score_design3_nf[[6]]),max(score_design4_nf[[6]]),
                         max(score_design5_nf[[6]]))
sizenf_score_after <- c(max(score_design1_nf[[3]]),max(score_design2_nf[[3]]),
                        max(score_design3_nf[[3]]),max(score_design4_nf[[3]]),
                        max(score_design5_nf[[3]]))
sizenf_score_random <- c(score_design1_nf_xi[2], score_design2_nf_xi[2],
                         score_design3_nf_xi[2], score_design4_nf_xi[2],
                         score_design5_nf_xi[2])
# Jeffreys
sizenf_Jbayes_before <- c(max(Jbayes_design1_nf[[6]]),max(Jbayes_design2_nf[[6]]),
                          max(Jbayes_design3_nf[[6]]),max(Jbayes_design4_nf[[6]]),
                          max(Jbayes_design5_nf[[6]]))
sizenf_Jbayes_after <- c(max(Jbayes_design1_nf[[3]]),max(Jbayes_design2_nf[[3]]),
                         max(Jbayes_design3_nf[[3]]),max(Jbayes_design4_nf[[3]]),
                         max(Jbayes_design5_nf[[3]]))
sizenf_Jbayes_random <- c(Jbayes_design1_nf_xi[2], Jbayes_design2_nf_xi[2],
                          Jbayes_design3_nf_xi[2], Jbayes_design4_nf_xi[2],
                          Jbayes_design5_nf_xi[2])
# Uniform
sizenf_Ubayes_before <- c(max(Ubayes_design1_nf[[6]]),max(Ubayes_design2_nf[[6]]),
                          max(Ubayes_design3_nf[[6]]),max(Ubayes_design4_nf[[6]]),
                          max(Ubayes_design5_nf[[6]]))
sizenf_Ubayes_after <- c(max(Ubayes_design1_nf[[3]]),max(Ubayes_design2_nf[[3]]),
                         max(Ubayes_design3_nf[[3]]),max(Ubayes_design4_nf[[3]]),
                         max(Ubayes_design5_nf[[3]]))
sizenf_Ubayes_random <- c(Ubayes_design1_nf_xi[2], Ubayes_design2_nf_xi[2],
                          Ubayes_design3_nf_xi[2], Ubayes_design4_nf_xi[2],
                          Ubayes_design5_nf_xi[2])
# logit normal 1
sizenf_Iln1_bayes_before <- c(max(Iln1_bayes_design1_nf[[6]]),max(Iln1_bayes_design2_nf[[6]]),
                              max(Iln1_bayes_design3_nf[[6]]),max(Iln1_bayes_design4_nf[[6]]),
                              max(Iln1_bayes_design5_nf[[6]]))
sizenf_Iln1_bayes_after <- c(max(Iln1_bayes_design1_nf[[3]]),max(Iln1_bayes_design2_nf[[3]]),
                             max(Iln1_bayes_design3_nf[[3]]),max(Iln1_bayes_design4_nf[[3]]),
                             max(Iln1_bayes_design5_nf[[3]]))
sizenf_Iln1_bayes_random <- c(Iln1_bayes_design1_nf_xi[2], Iln1_bayes_design2_nf_xi[2],
                              Iln1_bayes_design3_nf_xi[2], Iln1_bayes_design4_nf_xi[2],
                              Iln1_bayes_design5_nf_xi[2])
# logit normal 2
sizenf_Iln2_bayes_before <- c(max(Iln2_bayes_design1_nf[[6]]),max(Iln2_bayes_design2_nf[[6]]),
                              max(Iln2_bayes_design3_nf[[6]]),max(Iln2_bayes_design4_nf[[6]]),
                              max(Iln2_bayes_design5_nf[[6]]))
sizenf_Iln2_bayes_after <- c(max(Iln2_bayes_design1_nf[[3]]),max(Iln2_bayes_design2_nf[[3]]),
                             max(Iln2_bayes_design3_nf[[3]]),max(Iln2_bayes_design4_nf[[3]]),
                             max(Iln2_bayes_design5_nf[[3]]))
sizenf_Iln2_bayes_random <- c(Iln2_bayes_design1_nf_xi[2], Iln2_bayes_design2_nf_xi[2],
                              Iln2_bayes_design3_nf_xi[2], Iln2_bayes_design4_nf_xi[2],
                              Iln2_bayes_design5_nf_xi[2])
# corr logit normal
sizenf_Dln_bayes_before <- c(max(Dln_bayes_design1_nf[[6]]),max(Dln_bayes_design2_nf[[6]]),
                             max(Dln_bayes_design3_nf[[6]]),max(Dln_bayes_design4_nf[[6]]),
                             max(Dln_bayes_design5_nf[[6]]))
sizenf_Dln_bayes_after <- c(max(Dln_bayes_design1_nf[[3]]),max(Dln_bayes_design2_nf[[3]]),
                            max(Dln_bayes_design3_nf[[3]]),max(Dln_bayes_design4_nf[[3]]),
                            max(Dln_bayes_design5_nf[[3]]))
sizenf_Dln_bayes_random <- c(Dln_bayes_design1_nf_xi[2], Dln_bayes_design2_nf_xi[2],
                             Dln_bayes_design3_nf_xi[2], Dln_bayes_design4_nf_xi[2],
                             Dln_bayes_design5_nf_xi[2])
### Power
# unpooled z
powernf_zunpool_before <- c(power_calculate(zunpool_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                            power_calculate(zunpool_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                            power_calculate(zunpool_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                            power_calculate(zunpool_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                            power_calculate(zunpool_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_zunpool_after <- c(power_calculate(zunpool_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(zunpool_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(zunpool_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(zunpool_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(zunpool_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_zunpool_after_p2 <- c(zunpool_design1_nf_xi[1]*powerf_2nd(zunpool_design1_nf[[2]],zunpool_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                              zunpool_design2_nf_xi[1]*powerf_2nd(zunpool_design2_nf[[2]],zunpool_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                              zunpool_design3_nf_xi[1]*powerf_2nd(zunpool_design3_nf[[2]],zunpool_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                              zunpool_design4_nf_xi[1]*powerf_2nd(zunpool_design4_nf[[2]],zunpool_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                              zunpool_design5_nf_xi[1]*powerf_2nd(zunpool_design5_nf[[2]],zunpool_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_zunpool_random <- powernf_zunpool_after + powernf_zunpool_after_p2
# score
powernf_score_before <- c(power_calculate(score_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(score_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(score_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(score_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(score_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_score_after <- c(power_calculate(score_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                         power_calculate(score_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                         power_calculate(score_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                         power_calculate(score_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                         power_calculate(score_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_score_after_p2 <- c(score_design1_nf_xi[1]*powerf_2nd(score_design1_nf[[2]],score_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                            score_design2_nf_xi[1]*powerf_2nd(score_design2_nf[[2]],score_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                            score_design3_nf_xi[1]*powerf_2nd(score_design3_nf[[2]],score_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                            score_design4_nf_xi[1]*powerf_2nd(score_design4_nf[[2]],score_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                            score_design5_nf_xi[1]*powerf_2nd(score_design5_nf[[2]],score_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_score_random <- powernf_score_after + powernf_score_after_p2
# Jeffreys
powernf_Jbayes_before <- c(power_calculate(Jbayes_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Jbayes_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Jbayes_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Jbayes_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Jbayes_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Jbayes_after <- c(power_calculate(Jbayes_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Jbayes_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Jbayes_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Jbayes_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Jbayes_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Jbayes_after_p2 <- c(Jbayes_design1_nf_xi[1]*powerf_2nd(Jbayes_design1_nf[[2]],Jbayes_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                             Jbayes_design2_nf_xi[1]*powerf_2nd(Jbayes_design2_nf[[2]],Jbayes_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                             Jbayes_design3_nf_xi[1]*powerf_2nd(Jbayes_design3_nf[[2]],Jbayes_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                             Jbayes_design4_nf_xi[1]*powerf_2nd(Jbayes_design4_nf[[2]],Jbayes_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                             Jbayes_design5_nf_xi[1]*powerf_2nd(Jbayes_design5_nf[[2]],Jbayes_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_Jbayes_random <- powernf_Jbayes_after + powernf_Jbayes_after_p2
# Uniform
powernf_Ubayes_before <- c(power_calculate(Ubayes_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Ubayes_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Ubayes_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Ubayes_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                           power_calculate(Ubayes_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Ubayes_after <- c(power_calculate(Ubayes_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Ubayes_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Ubayes_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Ubayes_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                          power_calculate(Ubayes_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Ubayes_after_p2 <- c(Ubayes_design1_nf_xi[1]*powerf_2nd(Ubayes_design1_nf[[2]],Ubayes_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                             Ubayes_design2_nf_xi[1]*powerf_2nd(Ubayes_design2_nf[[2]],Ubayes_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                             Ubayes_design3_nf_xi[1]*powerf_2nd(Ubayes_design3_nf[[2]],Ubayes_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                             Ubayes_design4_nf_xi[1]*powerf_2nd(Ubayes_design4_nf[[2]],Ubayes_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                             Ubayes_design5_nf_xi[1]*powerf_2nd(Ubayes_design5_nf[[2]],Ubayes_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_Ubayes_random <- powernf_Ubayes_after + powernf_Ubayes_after_p2
# logit normal 1
powernf_Iln1_bayes_before <- c(power_calculate(Iln1_bayes_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln1_bayes_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln1_bayes_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln1_bayes_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln1_bayes_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Iln1_bayes_after <- c(power_calculate(Iln1_bayes_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln1_bayes_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln1_bayes_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln1_bayes_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln1_bayes_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Iln1_bayes_after_p2 <- c(Iln1_bayes_design1_nf_xi[1]*powerf_2nd(Iln1_bayes_design1_nf[[2]],Iln1_bayes_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                                 Iln1_bayes_design2_nf_xi[1]*powerf_2nd(Iln1_bayes_design2_nf[[2]],Iln1_bayes_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                                 Iln1_bayes_design3_nf_xi[1]*powerf_2nd(Iln1_bayes_design3_nf[[2]],Iln1_bayes_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                                 Iln1_bayes_design4_nf_xi[1]*powerf_2nd(Iln1_bayes_design4_nf[[2]],Iln1_bayes_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                                 Iln1_bayes_design5_nf_xi[1]*powerf_2nd(Iln1_bayes_design5_nf[[2]],Iln1_bayes_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_Iln1_bayes_random <- powernf_Iln1_bayes_after + powernf_Iln1_bayes_after_p2
# logit normal 2
powernf_Iln2_bayes_before <- c(power_calculate(Iln2_bayes_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln2_bayes_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln2_bayes_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln2_bayes_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                               power_calculate(Iln2_bayes_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Iln2_bayes_after <- c(power_calculate(Iln2_bayes_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln2_bayes_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln2_bayes_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln2_bayes_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Iln2_bayes_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Iln2_bayes_after_p2 <- c(Iln2_bayes_design1_nf_xi[1]*powerf_2nd(Iln2_bayes_design1_nf[[2]],Iln2_bayes_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                                 Iln2_bayes_design2_nf_xi[1]*powerf_2nd(Iln2_bayes_design2_nf[[2]],Iln2_bayes_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                                 Iln2_bayes_design3_nf_xi[1]*powerf_2nd(Iln2_bayes_design3_nf[[2]],Iln2_bayes_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                                 Iln2_bayes_design4_nf_xi[1]*powerf_2nd(Iln2_bayes_design4_nf[[2]],Iln2_bayes_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                                 Iln2_bayes_design5_nf_xi[1]*powerf_2nd(Iln2_bayes_design5_nf[[2]],Iln2_bayes_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_Iln2_bayes_random <- powernf_Iln2_bayes_after + powernf_Iln2_bayes_after_p2
# corr logit normal
powernf_Dln_bayes_before <- c(power_calculate(Dln_bayes_design1_nf[[5]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Dln_bayes_design2_nf[[5]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Dln_bayes_design3_nf[[5]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Dln_bayes_design4_nf[[5]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                              power_calculate(Dln_bayes_design5_nf[[5]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Dln_bayes_after <- c(power_calculate(Dln_bayes_design1_nf[[2]],nf_design[1]/4,nf_design[1]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                             power_calculate(Dln_bayes_design2_nf[[2]],nf_design[2]/4,nf_design[2]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                             power_calculate(Dln_bayes_design3_nf[[2]],nf_design[3]/4,nf_design[3]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                             power_calculate(Dln_bayes_design4_nf[[2]],nf_design[4]/4,nf_design[4]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06),
                             power_calculate(Dln_bayes_design5_nf[[2]],nf_design[5]/4,nf_design[5]/4*3,4.1/100,theta1 = 0.06, theta2 = 0.06))
powernf_Dln_bayes_after_p2 <- c(Dln_bayes_design1_nf_xi[1]*powerf_2nd(Dln_bayes_design1_nf[[2]],Dln_bayes_design1_nf_rr_part2,nf_design[1]/4,nf_design[1]/4*3,0,0.06),
                                Dln_bayes_design2_nf_xi[1]*powerf_2nd(Dln_bayes_design2_nf[[2]],Dln_bayes_design2_nf_rr_part2,nf_design[2]/4,nf_design[2]/4*3,0,0.06),
                                Dln_bayes_design3_nf_xi[1]*powerf_2nd(Dln_bayes_design3_nf[[2]],Dln_bayes_design3_nf_rr_part2,nf_design[3]/4,nf_design[3]/4*3,0,0.06),
                                Dln_bayes_design4_nf_xi[1]*powerf_2nd(Dln_bayes_design4_nf[[2]],Dln_bayes_design4_nf_rr_part2,nf_design[4]/4,nf_design[4]/4*3,0,0.06),
                                Dln_bayes_design5_nf_xi[1]*powerf_2nd(Dln_bayes_design5_nf[[2]],Dln_bayes_design5_nf_rr_part2,nf_design[5]/4,nf_design[5]/4*3,0,0.06))
powernf_Dln_bayes_random <- powernf_Dln_bayes_after + powernf_Dln_bayes_after_p2



# combined size
sizenf_before <- rbind(sizenf_zunpool_before, sizenf_score_before,
                       sizenf_Jbayes_before, sizenf_Ubayes_before, 
                       sizenf_Iln1_bayes_before, sizenf_Iln2_bayes_before, sizenf_Dln_bayes_before)
sizenf_after <- rbind(sizenf_zunpool_after, sizenf_score_after,
                      sizenf_Jbayes_after, sizenf_Ubayes_after, 
                      sizenf_Iln1_bayes_after, sizenf_Iln2_bayes_after, sizenf_Dln_bayes_after)
sizenf_random <- rbind(sizenf_zunpool_random, sizenf_score_random,
                       sizenf_Jbayes_random, sizenf_Ubayes_random, 
                       sizenf_Iln1_bayes_random, sizenf_Iln2_bayes_random, sizenf_Dln_bayes_random)
# combined power
power_before <- rbind(powernf_zunpool_before, powernf_score_before,
                      powernf_Jbayes_before, powernf_Ubayes_before, 
                      powernf_Iln1_bayes_before, powernf_Iln2_bayes_before, powernf_Dln_bayes_before)
power_after <- rbind(powernf_zunpool_after, powernf_score_after,
                     powernf_Jbayes_after, powernf_Ubayes_after,
                     powernf_Iln1_bayes_after, powernf_Iln2_bayes_after, powernf_Dln_bayes_after)
power_random <- rbind(powernf_zunpool_random, powernf_score_random,
                      powernf_Jbayes_random, powernf_Ubayes_random, 
                      powernf_Iln1_bayes_random, powernf_Iln2_bayes_random, powernf_Dln_bayes_random) # Note: this part is not plotted



