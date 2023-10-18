###########################################
##### Design of the Superiority Trial #####
###########################################

# Intended power: 0.9 at 0.15, 0.01; Allocation ratio 2:1 
sup_design <- seq(120, 180, 15)

##### Pooled Z #####
zpool_sup_time_overall <- c()
zpool_sup_time_overall <- rbind(zpool_sup_time_overall,
                                system.time(zpool_design1 <- 
                                              calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method="Pooled Z",gamma_initial=qnorm(0.975),
                                                                    theta=seq(0.001,0.999,0.001),alpha=0.025)))
zpool_sup_time_overall <- rbind(zpool_sup_time_overall,
                                system.time(zpool_design2 <- 
                                              calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method="Pooled Z",gamma_initial=qnorm(0.975),
                                                                    theta=seq(0.001,0.999,0.001),alpha=0.025)))
zpool_sup_time_overall <- rbind(zpool_sup_time_overall,
                                system.time(zpool_design3 <- 
                                              calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method="Pooled Z",gamma_initial=qnorm(0.975),
                                                                    theta=seq(0.001,0.999,0.001),alpha=0.025)))
zpool_sup_time_overall <- rbind(zpool_sup_time_overall,
                                system.time(zpool_design4 <- 
                                              calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method="Pooled Z",gamma_initial=qnorm(0.975),
                                                                    theta=seq(0.001,0.999,0.001),alpha=0.025)))
zpool_sup_time_overall <- rbind(zpool_sup_time_overall,
                                system.time(zpool_design5 <- 
                                              calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method="Pooled Z",gamma_initial=qnorm(0.975),
                                                                    theta=seq(0.001,0.999,0.001),alpha=0.025)))
# randomized part
zpool_design1_rr_part2 <- secondp_rd_freq(zpool_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Pooled Z")
zpool_design1_tIe_part2 <- powerf_2nd(zpool_design1[[2]], zpool_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
zpool_design1_xi <- binarys_xi(zpool_design1[[3]], zpool_design1_tIe_part2, 0.025, precision = 0.00001)

zpool_design2_rr_part2 <- secondp_rd_freq(zpool_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Pooled Z")
zpool_design2_tIe_part2 <- powerf_2nd(zpool_design2[[2]], zpool_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
zpool_design2_xi <- binarys_xi(zpool_design2[[3]], zpool_design2_tIe_part2, 0.025, precision = 0.00001)

zpool_design3_rr_part2 <- secondp_rd_freq(zpool_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Pooled Z")
zpool_design3_tIe_part2 <- powerf_2nd(zpool_design3[[2]], zpool_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
zpool_design3_xi <- binarys_xi(zpool_design3[[3]], zpool_design3_tIe_part2, 0.025, precision = 0.00001)

zpool_design4_rr_part2 <- secondp_rd_freq(zpool_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Pooled Z")
zpool_design4_tIe_part2 <- powerf_2nd(zpool_design4[[2]], zpool_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
zpool_design4_xi <- binarys_xi(zpool_design4[[3]], zpool_design4_tIe_part2, 0.025, precision = 0.00001)

zpool_design5_rr_part2 <- secondp_rd_freq(zpool_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Pooled Z")
zpool_design5_tIe_part2 <- powerf_2nd(zpool_design5[[2]], zpool_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
zpool_design5_xi <- binarys_xi(zpool_design5[[3]], zpool_design5_tIe_part2, 0.025, precision = 0.00001)



##### Unpooled Z #####
zunpool_sup_time_overall <- c()
zunpool_sup_time_overall <- rbind(zunpool_sup_time_overall,
                                  system.time(zunpool_design1 <- 
                                                calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method="Unpooled Z",gamma_initial=qnorm(0.975),
                                                                      theta=seq(0.001,0.999,0.001),alpha=0.025)))
zunpool_sup_time_overall <- rbind(zunpool_sup_time_overall,
                                  system.time(zunpool_design2 <- 
                                                calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method="Unpooled Z",gamma_initial=qnorm(0.975),
                                                                      theta=seq(0.001,0.999,0.001),alpha=0.025)))
zunpool_sup_time_overall <- rbind(zunpool_sup_time_overall,
                                  system.time(zunpool_design3 <- 
                                                calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method="Unpooled Z",gamma_initial=qnorm(0.975),
                                                                      theta=seq(0.001,0.999,0.001),alpha=0.025)))
zunpool_sup_time_overall <- rbind(zunpool_sup_time_overall,
                                  system.time(zunpool_design4 <- 
                                                calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method="Unpooled Z",gamma_initial=qnorm(0.975),
                                                                      theta=seq(0.001,0.999,0.001),alpha=0.025)))
zunpool_sup_time_overall <- rbind(zunpool_sup_time_overall,
                                  system.time(zunpool_design5 <- 
                                                calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method="Unpooled Z",gamma_initial=qnorm(0.975),
                                                                      theta=seq(0.001,0.999,0.001),alpha=0.025)))
# randomized part
zunpool_design1_rr_part2 <- secondp_rd_freq(zunpool_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Unpooled Z")
zunpool_design1_tIe_part2 <- powerf_2nd(zunpool_design1[[2]], zunpool_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
zunpool_design1_xi <- binarys_xi(zunpool_design1[[3]], zunpool_design1_tIe_part2, 0.025, precision = 0.00001)

zunpool_design2_rr_part2 <- secondp_rd_freq(zunpool_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Unpooled Z")
zunpool_design2_tIe_part2 <- powerf_2nd(zunpool_design2[[2]], zunpool_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
zunpool_design2_xi <- binarys_xi(zunpool_design2[[3]], zunpool_design2_tIe_part2, 0.025, precision = 0.00001)

zunpool_design3_rr_part2 <- secondp_rd_freq(zunpool_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Unpooled Z")
zunpool_design3_tIe_part2 <- powerf_2nd(zunpool_design3[[2]], zunpool_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
zunpool_design3_xi <- binarys_xi(zunpool_design3[[3]], zunpool_design3_tIe_part2, 0.025, precision = 0.00001)

zunpool_design4_rr_part2 <- secondp_rd_freq(zunpool_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Unpooled Z")
zunpool_design4_tIe_part2 <- powerf_2nd(zunpool_design4[[2]], zunpool_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
zunpool_design4_xi <- binarys_xi(zunpool_design4[[3]], zunpool_design4_tIe_part2, 0.025, precision = 0.00001)

zunpool_design5_rr_part2 <- secondp_rd_freq(zunpool_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Unpooled Z")
zunpool_design5_tIe_part2 <- powerf_2nd(zunpool_design5[[2]], zunpool_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
zunpool_design5_xi <- binarys_xi(zunpool_design5[[3]], zunpool_design5_tIe_part2, 0.025, precision = 0.00001)



##### Fisher's exact #####
FisherCp_sup_time_overall <- c()
FisherCp_sup_time_overall <- rbind(FisherCp_sup_time_overall,
                                   system.time(FisherCp_design1 <- 
                                                 calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method="Fisher's condition",gamma_initial=0.975,
                                                                       theta=seq(0.001,0.999,0.001),alpha=0.025)))
FisherCp_sup_time_overall <- rbind(FisherCp_sup_time_overall,
                                   system.time(FisherCp_design2 <- 
                                                 calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method="Fisher's condition",gamma_initial=0.975,
                                                                       theta=seq(0.001,0.999,0.001),alpha=0.025)))
FisherCp_sup_time_overall <- rbind(FisherCp_sup_time_overall,
                                   system.time(FisherCp_design3 <- 
                                                 calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method="Fisher's condition",gamma_initial=0.975,
                                                                       theta=seq(0.001,0.999,0.001),alpha=0.025)))
FisherCp_sup_time_overall <- rbind(FisherCp_sup_time_overall,
                                   system.time(FisherCp_design4 <- 
                                                 calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method="Fisher's condition",gamma_initial=0.975,
                                                                       theta=seq(0.001,0.999,0.001),alpha=0.025)))
FisherCp_sup_time_overall <- rbind(FisherCp_sup_time_overall,
                                   system.time(FisherCp_design5 <- 
                                                 calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method="Fisher's condition",gamma_initial=0.975,
                                                                       theta=seq(0.001,0.999,0.001),alpha=0.025)))
# randomized part
FisherCp_design1_rr_part2 <- secondp_rd_freq(FisherCp_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Fisher's condition")
FisherCp_design1_tIe_part2 <- powerf_2nd(FisherCp_design1[[2]], FisherCp_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
FisherCp_design1_xi <- binarys_xi(FisherCp_design1[[3]], FisherCp_design1_tIe_part2, 0.025, precision = 0.00001)

FisherCp_design2_rr_part2 <- secondp_rd_freq(FisherCp_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Fisher's condition")
FisherCp_design2_tIe_part2 <- powerf_2nd(FisherCp_design2[[2]], FisherCp_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
FisherCp_design2_xi <- binarys_xi(FisherCp_design2[[3]], FisherCp_design2_tIe_part2, 0.025, precision = 0.00001)

FisherCp_design3_rr_part2 <- secondp_rd_freq(FisherCp_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Fisher's condition")
FisherCp_design3_tIe_part2 <- powerf_2nd(FisherCp_design3[[2]], FisherCp_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
FisherCp_design3_xi <- binarys_xi(FisherCp_design3[[3]], FisherCp_design3_tIe_part2, 0.025, precision = 0.00001)

FisherCp_design4_rr_part2 <- secondp_rd_freq(FisherCp_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Fisher's condition")
FisherCp_design4_tIe_part2 <- powerf_2nd(FisherCp_design4[[2]], FisherCp_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
FisherCp_design4_xi <- binarys_xi(FisherCp_design4[[3]], FisherCp_design4_tIe_part2, 0.025, precision = 0.00001)

FisherCp_design5_rr_part2 <- secondp_rd_freq(FisherCp_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Fisher's condition")
FisherCp_design5_tIe_part2 <- powerf_2nd(FisherCp_design5[[2]], FisherCp_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
FisherCp_design5_xi <- binarys_xi(FisherCp_design5[[3]], FisherCp_design5_tIe_part2, 0.025, precision = 0.00001)



##### Bayesian (Jeffreys) #####
Jbayes_sup_time_overall <- c()
Jbayes_sup_time_overall <- rbind(Jbayes_sup_time_overall,
                                 system.time(Jbayes_design1 <- 
                                               calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_sup_time_overall <- rbind(Jbayes_sup_time_overall,
                                 system.time(Jbayes_design2 <- 
                                               calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_sup_time_overall <- rbind(Jbayes_sup_time_overall,
                                 system.time(Jbayes_design3 <- 
                                               calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_sup_time_overall <- rbind(Jbayes_sup_time_overall,
                                 system.time(Jbayes_design4 <- 
                                               calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0.5,0.5,0.5,0.5))))
Jbayes_sup_time_overall <- rbind(Jbayes_sup_time_overall,
                                 system.time(Jbayes_design5 <- 
                                               calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0.5,0.5,0.5,0.5))))
# randomized part
Jbayes_design1_rr_part2 <- secondp_rd_bayes(Jbayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design1_tIe_part2 <- powerf_2nd(Jbayes_design1[[2]], Jbayes_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
Jbayes_design1_xi <- binarys_xi(Jbayes_design1[[3]], Jbayes_design1_tIe_part2, 0.025, precision = 0.00001)

Jbayes_design2_rr_part2 <- secondp_rd_bayes(Jbayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design2_tIe_part2 <- powerf_2nd(Jbayes_design2[[2]], Jbayes_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
Jbayes_design2_xi <- binarys_xi(Jbayes_design2[[3]], Jbayes_design2_tIe_part2, 0.025, precision = 0.00001)

Jbayes_design3_rr_part2 <- secondp_rd_bayes(Jbayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design3_tIe_part2 <- powerf_2nd(Jbayes_design3[[2]], Jbayes_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
Jbayes_design3_xi <- binarys_xi(Jbayes_design3[[3]], Jbayes_design3_tIe_part2, 0.025, precision = 0.00001)

Jbayes_design4_rr_part2 <- secondp_rd_bayes(Jbayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design4_tIe_part2 <- powerf_2nd(Jbayes_design4[[2]], Jbayes_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
Jbayes_design4_xi <- binarys_xi(Jbayes_design4[[3]], Jbayes_design4_tIe_part2, 0.025, precision = 0.00001)

Jbayes_design5_rr_part2 <- secondp_rd_bayes(Jbayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Independent Beta", parameters=c(0.5,0.5,0.5,0.5))
Jbayes_design5_tIe_part2 <- powerf_2nd(Jbayes_design5[[2]], Jbayes_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
Jbayes_design5_xi <- binarys_xi(Jbayes_design5[[3]], Jbayes_design5_tIe_part2, 0.025, precision = 0.00001)



##### Bayes (Uniform) #####
Ubayes_sup_time_overall <- c()
Ubayes_sup_time_overall <- rbind(Ubayes_sup_time_overall,
                                 system.time(Ubayes_design1 <- 
                                               calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,1,1,1))))
Ubayes_sup_time_overall <- rbind(Ubayes_sup_time_overall,
                                 system.time(Ubayes_design2 <- 
                                               calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,1,1,1))))
Ubayes_sup_time_overall <- rbind(Ubayes_sup_time_overall,
                                 system.time(Ubayes_design3 <- 
                                               calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,1,1,1))))
Ubayes_sup_time_overall <- rbind(Ubayes_sup_time_overall,
                                 system.time(Ubayes_design4 <- 
                                               calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,1,1,1))))
Ubayes_sup_time_overall <- rbind(Ubayes_sup_time_overall,
                                 system.time(Ubayes_design5 <- 
                                               calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method=c("Bayes Posterior Probability", "Independent Beta"),gamma_initial=0.975,
                                                                     theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,1,1,1))))
# randomized part
Ubayes_design1_rr_part2 <- secondp_rd_bayes(Ubayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design1_tIe_part2 <- powerf_2nd(Ubayes_design1[[2]], Ubayes_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
Ubayes_design1_xi <- binarys_xi(Ubayes_design1[[3]], Ubayes_design1_tIe_part2, 0.025, precision = 0.00001)

Ubayes_design2_rr_part2 <- secondp_rd_bayes(Ubayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design2_tIe_part2 <- powerf_2nd(Ubayes_design2[[2]], Ubayes_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
Ubayes_design2_xi <- binarys_xi(Ubayes_design2[[3]], Ubayes_design2_tIe_part2, 0.025, precision = 0.00001)

Ubayes_design3_rr_part2 <- secondp_rd_bayes(Ubayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design3_tIe_part2 <- powerf_2nd(Ubayes_design3[[2]], Ubayes_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
Ubayes_design3_xi <- binarys_xi(Ubayes_design3[[3]], Ubayes_design3_tIe_part2, 0.025, precision = 0.00001)

Ubayes_design4_rr_part2 <- secondp_rd_bayes(Ubayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design4_tIe_part2 <- powerf_2nd(Ubayes_design4[[2]], Ubayes_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
Ubayes_design4_xi <- binarys_xi(Ubayes_design4[[3]], Ubayes_design4_tIe_part2, 0.025, precision = 0.00001)

Ubayes_design5_rr_part2 <- secondp_rd_bayes(Ubayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Independent Beta", parameters=c(1,1,1,1))
Ubayes_design5_tIe_part2 <- powerf_2nd(Ubayes_design5[[2]], Ubayes_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
Ubayes_design5_xi <- binarys_xi(Ubayes_design5[[3]], Ubayes_design5_tIe_part2, 0.025, precision = 0.00001)



##### Bayes (Logit N1) #####
Iln1_bayes_sup_time_overall <- c()
Iln1_bayes_sup_time_overall <- rbind(Iln1_bayes_sup_time_overall,
                                     system.time(Iln1_bayes_design1 <- 
                                                   calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,2,0,2,0))))
Iln1_bayes_sup_time_overall <- rbind(Iln1_bayes_sup_time_overall,
                                     system.time(Iln1_bayes_design2 <- 
                                                   calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,2,0,2,0))))
Iln1_bayes_sup_time_overall <- rbind(Iln1_bayes_sup_time_overall,
                                     system.time(Iln1_bayes_design3 <- 
                                                   calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,2,0,2,0))))
Iln1_bayes_sup_time_overall <- rbind(Iln1_bayes_sup_time_overall,
                                     system.time(Iln1_bayes_design4 <- 
                                                   calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,2,0,2,0))))
Iln1_bayes_sup_time_overall <- rbind(Iln1_bayes_sup_time_overall,
                                     system.time(Iln1_bayes_design5 <- 
                                                   calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,2,0,2,0))))
# randomized part
Iln1_bayes_design1_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design1_tIe_part2 <- powerf_2nd(Iln1_bayes_design1[[2]], Iln1_bayes_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
Iln1_bayes_design1_xi <- binarys_xi(Iln1_bayes_design1[[3]], Iln1_bayes_design1_tIe_part2, 0.025, precision = 0.00001)

Iln1_bayes_design2_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design2_tIe_part2 <- powerf_2nd(Iln1_bayes_design2[[2]], Iln1_bayes_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
Iln1_bayes_design2_xi <- binarys_xi(Iln1_bayes_design2[[3]], Iln1_bayes_design2_tIe_part2, 0.025, precision = 0.00001)

Iln1_bayes_design3_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design3_tIe_part2 <- powerf_2nd(Iln1_bayes_design3[[2]], Iln1_bayes_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
Iln1_bayes_design3_xi <- binarys_xi(Iln1_bayes_design3[[3]], Iln1_bayes_design3_tIe_part2, 0.025, precision = 0.00001)

Iln1_bayes_design4_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design4_tIe_part2 <- powerf_2nd(Iln1_bayes_design4[[2]], Iln1_bayes_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
Iln1_bayes_design4_xi <- binarys_xi(Iln1_bayes_design4[[3]], Iln1_bayes_design4_tIe_part2, 0.025, precision = 0.00001)

Iln1_bayes_design5_rr_part2 <- secondp_rd_bayes(Iln1_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Logit Normal",parameters=c(0,2,0,2,0))
Iln1_bayes_design5_tIe_part2 <- powerf_2nd(Iln1_bayes_design5[[2]], Iln1_bayes_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
Iln1_bayes_design5_xi <- binarys_xi(Iln1_bayes_design5[[3]], Iln1_bayes_design5_tIe_part2, 0.025, precision = 0.00001)



##### Bayes (Logit N2) #####
Iln2_bayes_sup_time_overall <- c()
Iln2_bayes_sup_time_overall <- rbind(Iln2_bayes_sup_time_overall,
                                     system.time(Iln2_bayes_design1 <- 
                                                   calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,3,0,3,0))))
Iln2_bayes_sup_time_overall <- rbind(Iln2_bayes_sup_time_overall,
                                     system.time(Iln2_bayes_design2 <- 
                                                   calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,3,0,3,0))))
Iln2_bayes_sup_time_overall <- rbind(Iln2_bayes_sup_time_overall,
                                     system.time(Iln2_bayes_design3 <- 
                                                   calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,3,0,3,0))))
Iln2_bayes_sup_time_overall <- rbind(Iln2_bayes_sup_time_overall,
                                     system.time(Iln2_bayes_design4 <- 
                                                   calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,3,0,3,0))))
Iln2_bayes_sup_time_overall <- rbind(Iln2_bayes_sup_time_overall,
                                     system.time(Iln2_bayes_design5 <- 
                                                   calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                         theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(0,3,0,3,0))))
# randomized part
Iln2_bayes_design1_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design1_tIe_part2 <- powerf_2nd(Iln2_bayes_design1[[2]], Iln2_bayes_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
Iln2_bayes_design1_xi <- binarys_xi(Iln2_bayes_design1[[3]], Iln2_bayes_design1_tIe_part2, 0.025, precision = 0.00001)

Iln2_bayes_design2_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design2_tIe_part2 <- powerf_2nd(Iln2_bayes_design2[[2]], Iln2_bayes_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
Iln2_bayes_design2_xi <- binarys_xi(Iln2_bayes_design2[[3]], Iln2_bayes_design2_tIe_part2, 0.025, precision = 0.00001)

Iln2_bayes_design3_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design3_tIe_part2 <- powerf_2nd(Iln2_bayes_design3[[2]], Iln2_bayes_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
Iln2_bayes_design3_xi <- binarys_xi(Iln2_bayes_design3[[3]], Iln2_bayes_design3_tIe_part2, 0.025, precision = 0.00001)

Iln2_bayes_design4_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design4_tIe_part2 <- powerf_2nd(Iln2_bayes_design4[[2]], Iln2_bayes_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
Iln2_bayes_design4_xi <- binarys_xi(Iln2_bayes_design4[[3]], Iln2_bayes_design4_tIe_part2, 0.025, precision = 0.00001)

Iln2_bayes_design5_rr_part2 <- secondp_rd_bayes(Iln2_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Logit Normal",parameters=c(0,3,0,3,0))
Iln2_bayes_design5_tIe_part2 <- powerf_2nd(Iln2_bayes_design5[[2]], Iln2_bayes_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
Iln2_bayes_design5_xi <- binarys_xi(Iln2_bayes_design5[[3]], Iln2_bayes_design5_tIe_part2, 0.025, precision = 0.00001)



##### Bayes (Logit NCorr) #####
Dln_bayes_sup_time_overall <- c()
Dln_bayes_sup_time_overall <- rbind(Dln_bayes_sup_time_overall,
                                    system.time(Dln_bayes_design1 <- 
                                                  calibration_algorithm(sup_design[1]/3*2,sup_design[1]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                        theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_sup_time_overall <- rbind(Dln_bayes_sup_time_overall,
                                    system.time(Dln_bayes_design2 <- 
                                                  calibration_algorithm(sup_design[2]/3*2,sup_design[2]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                        theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_sup_time_overall <- rbind(Dln_bayes_sup_time_overall,
                                    system.time(Dln_bayes_design3 <- 
                                                  calibration_algorithm(sup_design[3]/3*2,sup_design[3]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                        theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_sup_time_overall <- rbind(Dln_bayes_sup_time_overall,
                                    system.time(Dln_bayes_design4 <- 
                                                  calibration_algorithm(sup_design[4]/3*2,sup_design[4]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                        theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,3,-1,3,0.5))))
Dln_bayes_sup_time_overall <- rbind(Dln_bayes_sup_time_overall,
                                    system.time(Dln_bayes_design5 <- 
                                                  calibration_algorithm(sup_design[5]/3*2,sup_design[5]/3,0,method=c("Bayes Posterior Probability", "Logit Normal"),gamma_initial=0.975,
                                                                        theta=seq(0.001,0.999,0.001),alpha=0.025,parameters=c(1,3,-1,3,0.5))))
# randomized part
Dln_bayes_design1_rr_part2 <- secondp_rd_bayes(Dln_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design1_tIe_part2 <- powerf_2nd(Dln_bayes_design1[[2]], Dln_bayes_design1_rr_part2, sup_design[1]/3*2,sup_design[1]/3,0, seq(0.001,0.999,0.001))
Dln_bayes_design1_xi <- binarys_xi(Dln_bayes_design1[[3]], Dln_bayes_design1_tIe_part2, 0.025, precision = 0.00001)

Dln_bayes_design2_rr_part2 <- secondp_rd_bayes(Dln_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design2_tIe_part2 <- powerf_2nd(Dln_bayes_design2[[2]], Dln_bayes_design2_rr_part2, sup_design[2]/3*2,sup_design[2]/3,0, seq(0.001,0.999,0.001))
Dln_bayes_design2_xi <- binarys_xi(Dln_bayes_design2[[3]], Dln_bayes_design2_tIe_part2, 0.025, precision = 0.00001)

Dln_bayes_design3_rr_part2 <- secondp_rd_bayes(Dln_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design3_tIe_part2 <- powerf_2nd(Dln_bayes_design3[[2]], Dln_bayes_design3_rr_part2, sup_design[3]/3*2,sup_design[3]/3,0, seq(0.001,0.999,0.001))
Dln_bayes_design3_xi <- binarys_xi(Dln_bayes_design3[[3]], Dln_bayes_design3_tIe_part2, 0.025, precision = 0.00001)

Dln_bayes_design4_rr_part2 <- secondp_rd_bayes(Dln_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design4_tIe_part2 <- powerf_2nd(Dln_bayes_design4[[2]], Dln_bayes_design4_rr_part2, sup_design[4]/3*2,sup_design[4]/3,0, seq(0.001,0.999,0.001))
Dln_bayes_design4_xi <- binarys_xi(Dln_bayes_design4[[3]], Dln_bayes_design4_tIe_part2, 0.025, precision = 0.00001)

Dln_bayes_design5_rr_part2 <- secondp_rd_bayes(Dln_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,method="Logit Normal",parameters=c(1,3,-1,3,0.5))
Dln_bayes_design5_tIe_part2 <- powerf_2nd(Dln_bayes_design5[[2]], Dln_bayes_design5_rr_part2, sup_design[5]/3*2,sup_design[5]/3,0, seq(0.001,0.999,0.001))
Dln_bayes_design5_xi <- binarys_xi(Dln_bayes_design5[[3]], Dln_bayes_design5_tIe_part2, 0.025, precision = 0.00001)



##### Operating Characteristics #####
### Type I Error ###
# pooled z
sizesup_zpool_before <- c(max(zpool_design1[[6]]),max(zpool_design2[[6]]),
                          max(zpool_design3[[6]]),max(zpool_design4[[6]]),
                          max(zpool_design5[[6]]))
sizesup_zpool_after <- c(max(zpool_design1[[3]]),max(zpool_design2[[3]]),
                         max(zpool_design3[[3]]),max(zpool_design4[[3]]),
                         max(zpool_design5[[3]]))
sizesup_zpool_random <- c(zpool_design1_xi[2], zpool_design2_xi[2],
                          zpool_design3_xi[2], zpool_design4_xi[2],
                          zpool_design5_xi[2])
# unpooled z
sizesup_zunpool_before <- c(max(zunpool_design1[[6]]),max(zunpool_design2[[6]]),
                            max(zunpool_design3[[6]]),max(zunpool_design4[[6]]),
                            max(zunpool_design5[[6]]))
sizesup_zunpool_after <- c(max(zunpool_design1[[3]]),max(zunpool_design2[[3]]),
                           max(zunpool_design3[[3]]),max(zunpool_design4[[3]]),
                           max(zunpool_design5[[3]]))
sizesup_zunpool_random <- c(zunpool_design1_xi[2], zunpool_design2_xi[2],
                            zunpool_design3_xi[2], zunpool_design4_xi[2],
                            zunpool_design5_xi[2])
# fisher
sizesup_FisherCp_before <- c(max(FisherCp_design1[[6]]),max(FisherCp_design2[[6]]),
                             max(FisherCp_design3[[6]]),max(FisherCp_design4[[6]]),
                             max(FisherCp_design5[[6]]))
sizesup_FisherCp_after <- c(max(FisherCp_design1[[3]]),max(FisherCp_design2[[3]]),
                            max(FisherCp_design3[[3]]),max(FisherCp_design4[[3]]),
                            max(FisherCp_design5[[3]]))
sizesup_FisherCp_random <- c(FisherCp_design1_xi[2], FisherCp_design2_xi[2],
                             FisherCp_design3_xi[2], FisherCp_design4_xi[2],
                             FisherCp_design5_xi[2])
# Jeffreys
sizesup_Jbayes_before <- c(max(Jbayes_design1[[6]]),max(Jbayes_design2[[6]]),
                           max(Jbayes_design3[[6]]),max(Jbayes_design4[[6]]),
                           max(Jbayes_design5[[6]]))
sizesup_Jbayes_after <- c(max(Jbayes_design1[[3]]),max(Jbayes_design2[[3]]),
                          max(Jbayes_design3[[3]]),max(Jbayes_design4[[3]]),
                          max(Jbayes_design5[[3]]))
sizesup_Jbayes_random <- c(Jbayes_design1_xi[2], Jbayes_design2_xi[2],
                           Jbayes_design3_xi[2], Jbayes_design4_xi[2],
                           Jbayes_design5_xi[2])
# Uniform
sizesup_Ubayes_before <- c(max(Ubayes_design1[[6]]),max(Ubayes_design2[[6]]),
                           max(Ubayes_design3[[6]]),max(Ubayes_design4[[6]]),
                           max(Ubayes_design5[[6]]))
sizesup_Ubayes_after <- c(max(Ubayes_design1[[3]]),max(Ubayes_design2[[3]]),
                          max(Ubayes_design3[[3]]),max(Ubayes_design4[[3]]),
                          max(Ubayes_design5[[3]]))
sizesup_Ubayes_random <- c(Ubayes_design1_xi[2], Ubayes_design2_xi[2],
                           Ubayes_design3_xi[2], Ubayes_design4_xi[2],
                           Ubayes_design5_xi[2])
# logit normal 1
sizesup_Iln1_bayes_before <- c(max(Iln1_bayes_design1[[6]]),max(Iln1_bayes_design2[[6]]),
                               max(Iln1_bayes_design3[[6]]),max(Iln1_bayes_design4[[6]]),
                               max(Iln1_bayes_design5[[6]]))
sizesup_Iln1_bayes_after <- c(max(Iln1_bayes_design1[[3]]),max(Iln1_bayes_design2[[3]]),
                              max(Iln1_bayes_design3[[3]]),max(Iln1_bayes_design4[[3]]),
                              max(Iln1_bayes_design5[[3]]))
sizesup_Iln1_bayes_random <- c(Iln1_bayes_design1_xi[2], Iln1_bayes_design2_xi[2],
                               Iln1_bayes_design3_xi[2], Iln1_bayes_design4_xi[2],
                               Iln1_bayes_design5_xi[2])
# logit normal 2
sizesup_Iln2_bayes_before <- c(max(Iln2_bayes_design1[[6]]),max(Iln2_bayes_design2[[6]]),
                               max(Iln2_bayes_design3[[6]]),max(Iln2_bayes_design4[[6]]),
                               max(Iln2_bayes_design5[[6]]))
sizesup_Iln2_bayes_after <- c(max(Iln2_bayes_design1[[3]]),max(Iln2_bayes_design2[[3]]),
                              max(Iln2_bayes_design3[[3]]),max(Iln2_bayes_design4[[3]]),
                              max(Iln2_bayes_design5[[3]]))
sizesup_Iln2_bayes_random <- c(Iln2_bayes_design1_xi[2], Iln2_bayes_design2_xi[2],
                               Iln2_bayes_design3_xi[2], Iln2_bayes_design4_xi[2],
                               Iln2_bayes_design5_xi[2])
# corr logit normal
sizesup_Dln_bayes_before <- c(max(Dln_bayes_design1[[6]]),max(Dln_bayes_design2[[6]]),
                              max(Dln_bayes_design3[[6]]),max(Dln_bayes_design4[[6]]),
                              max(Dln_bayes_design5[[6]]))
sizesup_Dln_bayes_after <- c(max(Dln_bayes_design1[[3]]),max(Dln_bayes_design2[[3]]),
                             max(Dln_bayes_design3[[3]]),max(Dln_bayes_design4[[3]]),
                             max(Dln_bayes_design5[[3]]))
sizesup_Dln_bayes_random <- c(Dln_bayes_design1_xi[2], Dln_bayes_design2_xi[2],
                              Dln_bayes_design3_xi[2], Dln_bayes_design4_xi[2],
                              Dln_bayes_design5_xi[2])
### Power ###
# pooled z
powersup_zpool_before <- c(power_calculate(zpool_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(zpool_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(zpool_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(zpool_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(zpool_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_zpool_after <- c(power_calculate(zpool_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                          power_calculate(zpool_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                          power_calculate(zpool_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                          power_calculate(zpool_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                          power_calculate(zpool_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_zpool_after_p2 <- c(zpool_design1_xi[1]*powerf_2nd(zpool_design1[[2]],zpool_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                             zpool_design2_xi[1]*powerf_2nd(zpool_design2[[2]],zpool_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                             zpool_design3_xi[1]*powerf_2nd(zpool_design3[[2]],zpool_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                             zpool_design4_xi[1]*powerf_2nd(zpool_design4[[2]],zpool_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                             zpool_design5_xi[1]*powerf_2nd(zpool_design5[[2]],zpool_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_zpool_random <- powersup_zpool_after + powersup_zpool_after_p2
# unpooled z
powersup_zunpool_before <- c(power_calculate(zunpool_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(zunpool_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(zunpool_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(zunpool_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(zunpool_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_zunpool_after <- c(power_calculate(zunpool_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(zunpool_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(zunpool_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(zunpool_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(zunpool_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_zunpool_after_p2 <- c(zunpool_design1_xi[1]*powerf_2nd(zunpool_design1[[2]],zunpool_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                               zunpool_design2_xi[1]*powerf_2nd(zunpool_design2[[2]],zunpool_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                               zunpool_design3_xi[1]*powerf_2nd(zunpool_design3[[2]],zunpool_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                               zunpool_design4_xi[1]*powerf_2nd(zunpool_design4[[2]],zunpool_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                               zunpool_design5_xi[1]*powerf_2nd(zunpool_design5[[2]],zunpool_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_zunpool_random <- powersup_zunpool_after + powersup_zunpool_after_p2
# fisher
powersup_FisherCp_before <- c(power_calculate(FisherCp_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(FisherCp_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(FisherCp_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(FisherCp_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(FisherCp_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_FisherCp_after <- c(power_calculate(FisherCp_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(FisherCp_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(FisherCp_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(FisherCp_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                             power_calculate(FisherCp_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_FisherCp_after_p2 <- c(FisherCp_design1_xi[1]*powerf_2nd(FisherCp_design1[[2]],FisherCp_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                                FisherCp_design2_xi[1]*powerf_2nd(FisherCp_design2[[2]],FisherCp_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                                FisherCp_design3_xi[1]*powerf_2nd(FisherCp_design3[[2]],FisherCp_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                                FisherCp_design4_xi[1]*powerf_2nd(FisherCp_design4[[2]],FisherCp_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                                FisherCp_design5_xi[1]*powerf_2nd(FisherCp_design5[[2]],FisherCp_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_FisherCp_random <- powersup_FisherCp_after + powersup_FisherCp_after_p2
# Jeffreys
powersup_Jbayes_before <- c(power_calculate(Jbayes_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Jbayes_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Jbayes_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Jbayes_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Jbayes_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Jbayes_after <- c(power_calculate(Jbayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Jbayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Jbayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Jbayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Jbayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Jbayes_after_p2 <- c(Jbayes_design1_xi[1]*powerf_2nd(Jbayes_design1[[2]],Jbayes_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                              Jbayes_design2_xi[1]*powerf_2nd(Jbayes_design2[[2]],Jbayes_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                              Jbayes_design3_xi[1]*powerf_2nd(Jbayes_design3[[2]],Jbayes_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                              Jbayes_design4_xi[1]*powerf_2nd(Jbayes_design4[[2]],Jbayes_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                              Jbayes_design5_xi[1]*powerf_2nd(Jbayes_design5[[2]],Jbayes_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_Jbayes_random <- powersup_Jbayes_after + powersup_Jbayes_after_p2
# Uniform
powersup_Ubayes_before <- c(power_calculate(Ubayes_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Ubayes_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Ubayes_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Ubayes_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                            power_calculate(Ubayes_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Ubayes_after <- c(power_calculate(Ubayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Ubayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Ubayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Ubayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                           power_calculate(Ubayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Ubayes_after_p2 <- c(Ubayes_design1_xi[1]*powerf_2nd(Ubayes_design1[[2]],Ubayes_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                              Ubayes_design2_xi[1]*powerf_2nd(Ubayes_design2[[2]],Ubayes_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                              Ubayes_design3_xi[1]*powerf_2nd(Ubayes_design3[[2]],Ubayes_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                              Ubayes_design4_xi[1]*powerf_2nd(Ubayes_design4[[2]],Ubayes_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                              Ubayes_design5_xi[1]*powerf_2nd(Ubayes_design5[[2]],Ubayes_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_Ubayes_random <- powersup_Ubayes_after + powersup_Ubayes_after_p2
# logit normal 1
powersup_Iln1_bayes_before <- c(power_calculate(Iln1_bayes_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln1_bayes_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln1_bayes_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln1_bayes_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln1_bayes_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Iln1_bayes_after <- c(power_calculate(Iln1_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln1_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln1_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln1_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln1_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Iln1_bayes_after_p2 <- c(Iln1_bayes_design1_xi[1]*powerf_2nd(Iln1_bayes_design1[[2]],Iln1_bayes_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                                  Iln1_bayes_design2_xi[1]*powerf_2nd(Iln1_bayes_design2[[2]],Iln1_bayes_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                                  Iln1_bayes_design3_xi[1]*powerf_2nd(Iln1_bayes_design3[[2]],Iln1_bayes_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                                  Iln1_bayes_design4_xi[1]*powerf_2nd(Iln1_bayes_design4[[2]],Iln1_bayes_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                                  Iln1_bayes_design5_xi[1]*powerf_2nd(Iln1_bayes_design5[[2]],Iln1_bayes_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_Iln1_bayes_random <- powersup_Iln1_bayes_after + powersup_Iln1_bayes_after_p2
# logit normal 2
powersup_Iln2_bayes_before <- c(power_calculate(Iln2_bayes_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln2_bayes_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln2_bayes_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln2_bayes_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                                power_calculate(Iln2_bayes_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Iln2_bayes_after <- c(power_calculate(Iln2_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln2_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln2_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln2_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Iln2_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Iln2_bayes_after_p2 <- c(Iln2_bayes_design1_xi[1]*powerf_2nd(Iln2_bayes_design1[[2]],Iln2_bayes_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                                  Iln2_bayes_design2_xi[1]*powerf_2nd(Iln2_bayes_design2[[2]],Iln2_bayes_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                                  Iln2_bayes_design3_xi[1]*powerf_2nd(Iln2_bayes_design3[[2]],Iln2_bayes_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                                  Iln2_bayes_design4_xi[1]*powerf_2nd(Iln2_bayes_design4[[2]],Iln2_bayes_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                                  Iln2_bayes_design5_xi[1]*powerf_2nd(Iln2_bayes_design5[[2]],Iln2_bayes_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_Iln2_bayes_random <- powersup_Iln2_bayes_after + powersup_Iln2_bayes_after_p2
# corr logit normal
powersup_Dln_bayes_before <- c(power_calculate(Dln_bayes_design1[[5]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Dln_bayes_design2[[5]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Dln_bayes_design3[[5]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Dln_bayes_design4[[5]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                               power_calculate(Dln_bayes_design5[[5]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Dln_bayes_after <- c(power_calculate(Dln_bayes_design1[[2]],sup_design[1]/3*2,sup_design[1]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(Dln_bayes_design2[[2]],sup_design[2]/3*2,sup_design[2]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(Dln_bayes_design3[[2]],sup_design[3]/3*2,sup_design[3]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(Dln_bayes_design4[[2]],sup_design[4]/3*2,sup_design[4]/3,0,theta1 = 0.15, theta2 = 0.01),
                              power_calculate(Dln_bayes_design5[[2]],sup_design[5]/3*2,sup_design[5]/3,0,theta1 = 0.15, theta2 = 0.01))
powersup_Dln_bayes_after_p2 <- c(Dln_bayes_design1_xi[1]*powerf_2nd(Dln_bayes_design1[[2]],Dln_bayes_design1_rr_part2,sup_design[1]/3*2,sup_design[1]/3,-0.14,0.15),
                                 Dln_bayes_design2_xi[1]*powerf_2nd(Dln_bayes_design2[[2]],Dln_bayes_design2_rr_part2,sup_design[2]/3*2,sup_design[2]/3,-0.14,0.15),
                                 Dln_bayes_design3_xi[1]*powerf_2nd(Dln_bayes_design3[[2]],Dln_bayes_design3_rr_part2,sup_design[3]/3*2,sup_design[3]/3,-0.14,0.15),
                                 Dln_bayes_design4_xi[1]*powerf_2nd(Dln_bayes_design4[[2]],Dln_bayes_design4_rr_part2,sup_design[4]/3*2,sup_design[4]/3,-0.14,0.15),
                                 Dln_bayes_design5_xi[1]*powerf_2nd(Dln_bayes_design5[[2]],Dln_bayes_design5_rr_part2,sup_design[5]/3*2,sup_design[5]/3,-0.14,0.15))
powersup_Dln_bayes_random <- powersup_Dln_bayes_after + powersup_Dln_bayes_after_p2



# combined size
sizesup_before <- rbind(sizesup_zpool_before, sizesup_zunpool_before, sizesup_FisherCp_before,
                        sizesup_Jbayes_before, sizesup_Ubayes_before, sizesup_Iln1_bayes_before,
                        sizesup_Iln2_bayes_before, sizesup_Dln_bayes_before)
sizesup_after <- rbind(sizesup_zpool_after, sizesup_zunpool_after, sizesup_FisherCp_after,
                       sizesup_Jbayes_after, sizesup_Ubayes_after, sizesup_Iln1_bayes_after,
                       sizesup_Iln2_bayes_after, sizesup_Dln_bayes_after)
sizesup_random <- rbind(sizesup_zpool_random, sizesup_zunpool_random, sizesup_FisherCp_random,
                        sizesup_Jbayes_random, sizesup_Ubayes_random, sizesup_Iln1_bayes_random,
                        sizesup_Iln2_bayes_random, sizesup_Dln_bayes_random)
# combined power
power_before <- rbind(powersup_zpool_before, powersup_zunpool_before, powersup_FisherCp_before,
                      powersup_Jbayes_before, powersup_Ubayes_before, powersup_Iln1_bayes_before,
                      powersup_Iln2_bayes_before, powersup_Dln_bayes_before)
power_after <- rbind(powersup_zpool_after, powersup_zunpool_after, powersup_FisherCp_after,
                     powersup_Jbayes_after, powersup_Ubayes_after, powersup_Iln1_bayes_after,
                     powersup_Iln2_bayes_after, powersup_Dln_bayes_after)
power_random <- rbind(powersup_zpool_random, powersup_zunpool_random, powersup_FisherCp_random,
                      powersup_Jbayes_random, powersup_Ubayes_random, powersup_Iln1_bayes_random,
                      powersup_Iln2_bayes_random, powersup_Dln_bayes_random) # Note: this part is not plotted


