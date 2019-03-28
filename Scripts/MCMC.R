# Reproducible code for "Bayesian atmospheric tomography for detection and quantification of methane emissions: Application to data from the 2015 Ginninderra release experiment​" by Cartwright et al.
# Copyright (c) 2019 Laura Cartwright
# Author: Laura Cartwright (lcartwri@uow.edu.au)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

mcmc_run <- function(mcmc_data, H, samples, iterates, burnin, thins){
  
  data_size <- nrow(mcmc_data)
  air_temp <- mcmc_data$air_temp
  air_press <- mcmc_data$air_pressure
  observeds <- mcmc_data$bg_conc
  WindSp <- mcmc_data$wind_speed
  
  # Create variable which equals wind speed if wind speed < 1, or 1 for wind speed >= 1
  if (tau_WS == 1){
    WindSp_hat <- c()
    WindSp_hat[mcmc_data$wind_speed < 1] <- mcmc_data$wind_speed[mcmc_data$wind_speed < 1]
    WindSp_hat[mcmc_data$wind_speed >= 1] <- 1
  } else {
    WindSp_hat <- rep(1, data_size)
  }
  
  # Create factor variable for combinations of variables for each tau
  if (Tau_track == 1){
    inst_grp <- as.factor(mcmc_data$inst_group) #careful with number of levels here
    if (stab_track == 1){
      second_group <- as.factor(mcmc_data$Stability)
    } else {
      second_group <- rep("1", data_size)
    }
    
    classification <- c()
    for (i in 1:data_size){
      grp <- inst_grp[i]
      dy <- second_group[i]
      classification[i] <- paste0(grp,":",dy)
    }
    classification <- as.factor(classification)
    igroups <- length(levels(classification))
  } else {
    classification <- as.factor(rep("one", data_size))
    igroups <- 1
  }
  
  ################################################################################
  
  Qsd <- 0.1 # initial SD for Q proposal
  SD <- rep(0.1, data_size) # initial SD for wind speed proposals
  tausd <- rep(0.1, igroups) # initial SD for tau proposals
  sigy_sd <- 1 # initial SD for sigy parameter proposals
  sigz_sd <- 1 # initial SD for sigz parameter proposals
  Q_prior_sd <- 1.5
  
  ################################################################################
  
  # Create predicted concentrations
  
  stability <- matrix(nrow = data_size, ncol = 4)
  for (i in 1:data_size){
    stability[i, ] <- stability_classes[ , mcmc_data$Stability[i]]
  }
  
  which_pts <- which(is.na(mcmc_data$x2))
  which_pth <- which(is.na(mcmc_data$x2) == FALSE)
  
  pts <- mcmc_data[which_pts, ]
  pth <- mcmc_data[which_pth, ]
  
  # Build pointwise line for path measurements
  ss <- seq(0, 1, by = 1 / (samples - 1))
  x_coords <- matrix(nrow = nrow(pth), ncol = samples)
  y_coords <- matrix(nrow = nrow(pth), ncol = samples)
  if (nrow(pth) >= 1){
      for (i in 1:nrow(pth)){
        x_coords[i, ] <- pth$x2[i] * ss + (1 - ss) * pth$x1[i]
        y_coords[i, ] <- pth$y2[i] * ss + (1 - ss) * pth$y1[i]
      }
  }
  
  # Calculate first set of predicted values & set initial values for sigma parameters
  if (sigy_track == 1){
    sigy <- c()
    sigy[1] <- 1 # initial value for sigy
    acceptance_sigy <- c()
    acceptance_sigy[1] <- 1
  } else {
    sigy <- rep(1, iterates)
    acceptance_sigy <- rep(1, iterates)
  }
  
  if (sigz_track == 1){
    sigz <- c()
    sigz[1] <- 1 # initial value for sigy
    acceptance_sigz <- c()
    acceptance_sigz[1] <- 1
  } else {
    sigz <- rep(1, iterates)
    acceptance_sigz <- rep(1, iterates)
  }
  
  predicteds <- c()
  
    pts_preds <- gaussian_plume_mat(as.matrix(pts$x1), as.matrix(pts$y1), pts$z, H, stability[which_pts, ], sigy[1], sigz[1])
    pth_preds <- gaussian_plume_mat(x_coords, y_coords, pth$z, H, stability[which_pth, ], sigy[1], sigz[1])
    predicteds[which_pts] <- pts_preds
    predicteds[which_pth] <- pth_preds
    
  
  ################################################################################
  
  # Set initial/start values
  
  # Q
  Q <- c()
  Q[1] <- 1
  acceptance <- c() # Keep track of acceptance rate of Q
  acceptance[1] <- 1
  
  # Taus
  if (Tau_track == 0){
    Tau <- 4
    tau <- matrix(nrow = iterates, ncol = igroups)
  } else {
    tau <- matrix(nrow = iterates, ncol = igroups)
    tau[1, ] <- rep(1, igroups)
    acceptance_tau <- matrix(nrow = 50, ncol = igroups)
  }
  
  
  ################################################################################
  
  # Run the MCMC

  for (i in 2:iterates){
    if (i == 2){
        print(Sys.time(), quote = FALSE) # print start time
    }
    if (i %% (iterates * 0.1) == 0){
      print(paste0((i / iterates) * 100, "% done"), quote = FALSE)
    }
    
    
    ## #######
    
    predicteds_tau <- (Q[i - 1] / WindSp) * predicteds
    predicteds_tau <- gm3toppm(predicteds_tau, air_temp, air_press, 16.04)
    
    half_diff_sq1 <- c() #observeds - predicteds squared and halved (needed only for normal likelihood i think)
    sizes <- c() # Number of obvs in each tau group
    
    for (j in 1:igroups){
        sizes[j] <- length(mcmc_data$bg_conc[classification == levels(classification)[j]])
      half_diff_sq1[j] <- sum((observeds[classification == levels(classification)[j]] - 
                                 predicteds_tau[classification == levels(classification)[j]]) ^ 2 / 2)
    }
    
    ##### Tau
    if (Tau_track == 1){
      if (liketype == "laplace"){ # Need to accept/reject if using Laplace likelihood
        if (i %% 10 == 0 & i >= 51 & i < 1000){
          tau_acc_r8 <- colMeans(acceptance_tau)
          tausd[tau_acc_r8 < 0.1] <- tausd[tau_acc_r8 < 0.1] * 0.9 # adapt tau SD on proposals
          tausd[tau_acc_r8 > 0.8] <- tausd[tau_acc_r8 > 0.8] * 1.1
        } else if (i %% 100 == 0 & i >= 1000 & i <= 20000){
          tau_acc_r8 <- colMeans(acceptance_tau)
          tausd[tau_acc_r8 < 0.1] <- tausd[tau_acc_r8 < 0.1] * 0.9
          tausd[tau_acc_r8 > 0.8] <- tausd[tau_acc_r8 > 0.8] * 1.1
        } else if (i == 20001) {
          acc_tau <- colMeans(acceptance_tau)
          print(paste0("# of tau acceptance rates outside (0.1, 0.8) at completion of adaptation: ", 
                       length(acc_tau[acc_tau < 0.1 | acc_tau > 0.8])), quote = FALSE)
        }
        
        # Sample
        tau_prop <- rnorm(igroups, mean = tau[(i - 1), ], sd = tausd) # propose values
        
        # Calculate log ratio
        log_prop_ratio_tau <- (log_tau(tau_prop, 1.058, 0.621, classification) + 
                                 log_like_tau(observeds, predicteds_tau, tau_prop,
                                              classification, sizes, tau_WS, WindSp_hat)) -
                              (log_tau(tau[(i - 1), ], 1.058, 0.621, classification) + 
                                 log_like_tau(observeds, predicteds_tau, tau[(i - 1), ],
                                              classification, sizes, tau_WS, WindSp_hat))
      
        # Exp(log ratio)
        prop_ratio_tau <- exp(log_prop_ratio_tau)
        
        # Accept/reject
        uu <- runif(igroups)
        accept_tau <- uu < prop_ratio_tau
        new_tau <- tau_prop * accept_tau + tau[(i - 1), ] * (!accept_tau)
        tau[i, ] <- new_tau
        acceptance_tau[1:49, ] <- acceptance_tau[2:50, ] 
        acceptance_tau[50, ] <- accept_tau
        
      } else if (liketype == "normal"){
        tau[i, ] <- rgamma(igroups, shape = 1.058 + (sizes / 2), rate = 0.621 + half_diff_sq1)
      }
    } else {
      tau[i, ] <- Tau
    }
    
    
    ######## new sigy
    if (sigy_track == 1){
      if ((i %% 40 == 0 & i >= 50 & i < 1000) | (i == 50)){
        if (mean(acceptance_sigy[(i - 40):(i - 1)]) < 0.1){ # Adapt sigy SD on proposal
          sigy_sd <- sigy_sd * 0.9
        } else if (mean(acceptance_sigy[(i - 40):(i - 1)]) > 0.8) {
          sigy_sd <- sigy_sd * 1.1
        }
      } else if (i %% 100 == 0 & i >= 1000 & i <= 20000){
        if (mean(acceptance_sigy[(i - 100):(i - 1)]) < 0.1){
          sigy_sd <- sigy_sd * 0.9
        } else if (mean(acceptance_sigy[(i - 100):(i - 1)]) > 0.8) {
          sigy_sd <- sigy_sd * 1.1
        }
        if (i == 20000){
          print(paste0("omega_y acceptance rate at completion of adaptation: ", 
                       mean(acceptance_sigy)), quote = FALSE)
        }
      }
      
      sigy_prop <- rnorm(1, sigy[i - 1], sigy_sd)
      
      # Create new predicteds with proposal
      
      predicteds1 <- c() # predicteds using proposed sigy
      
      pts_preds <- gaussian_plume_mat(as.matrix(pts$x1), as.matrix(pts$y1), pts$z, H, stability[which_pts, ], sigy_prop, sigz[i - 1])
      pth_preds <- gaussian_plume_mat(x_coords, y_coords, pth$z, H, stability[which_pth, ], sigy_prop, sigz[i - 1])
      predicteds1[which_pts] <- pts_preds
      predicteds1[which_pth] <- pth_preds
      
      # Create transformed predicteds
    
      predicteds_prop <- (Q[i - 1] / WindSp) * predicteds1
      predicteds_prop <- gm3toppm(predicteds_prop, air_temp, air_press, 16.04)
      
      # and old ones
      
      predicteds_old <- (Q[i - 1] / WindSp) * predicteds
      predicteds_old <- gm3toppm(predicteds_old, air_temp, air_press, 16.04)
      
      # Calculate log ratio
      
      log_prop_ratio_sigy <- (dgamma(sigy_prop, shape = 1.608356, rate = 0.7360556, log = TRUE) +
                              log_like(observeds, predicteds_prop, tau[i, ],
                                  data_size, classification, liketype, tau_WS, WindSp_hat)) -
                            (dgamma(sigy[i - 1], shape = 1.608356, rate = 0.7360556, log = TRUE) +
                              log_like(observeds, predicteds_old, tau[i, ],
                                  data_size, classification, liketype, tau_WS, WindSp_hat))
      
      # Exp(log ratio)
      prop_ratio_sigy <- exp(log_prop_ratio_sigy)
      
      # Accept/reject
      if (prop_ratio_sigy >= 1){
          sigy[i] <- sigy_prop
          predicteds <- predicteds1
          acceptance_sigy[i] <- 1
      } else {
          uu <- runif(1)
          if (uu <= prop_ratio_sigy){
              sigy[i] <- sigy_prop
              predicteds <- predicteds1
              acceptance_sigy[i] <- 1
          } else {
              sigy[i] <- sigy[i - 1]
              predicteds <- predicteds
              acceptance_sigy[i] <- 0
          }
      }
    }
    
    ######## new sigz
    
    if (sigz_track == 1){
      if ((i %% 40 == 0 & i >= 50 & i < 1000) | (i == 50)){
        if (mean(acceptance_sigz[(i - 40):(i - 1)]) < 0.1){ # Adapt sigz SD on proposal
          sigz_sd <- sigz_sd * 0.9
        } else if (mean(acceptance_sigz[(i - 40):(i - 1)]) > 0.8) {
          sigz_sd <- sigz_sd * 1.1
        }
      } else if (i %% 100 == 0 & i >= 1000 & i <= 20000){
        if (mean(acceptance_sigz[(i - 100):(i - 1)]) < 0.1){
          sigz_sd <- sigz_sd * 0.9
        } else if (mean(acceptance_sigz[(i - 100):(i - 1)]) > 0.8) {
          sigz_sd <- sigz_sd * 1.1
        }
        if (i == 20000){
          print(paste0("omega_z acceptance rate at completion of adaptation: ", 
                       mean(acceptance_sigz)), quote = FALSE)
        }
      }
      
      sigz_prop <- rnorm(1, sigz[i - 1], sigz_sd)
      
      # Create new predicteds with proposal
      
      predicteds1 <- c() # predicteds using proposed sigz
      
      pts_preds <- gaussian_plume_mat(as.matrix(pts$x1), as.matrix(pts$y1), pts$z, H, stability[which_pts, ], sigy[i], sigz_prop)
      pth_preds <- gaussian_plume_mat(x_coords, y_coords, pth$z, H, stability[which_pth, ], sigy[i], sigz_prop)
      predicteds1[which_pts] <- pts_preds
      predicteds1[which_pth] <- pth_preds
      
      # Create transformed predicteds
      
      predicteds_prop <- (Q[i - 1] / WindSp) * predicteds1
      predicteds_prop <- gm3toppm(predicteds_prop, air_temp, air_press, 16.04)
      
      # and old ones
      
      predicteds_old <- (Q[i - 1] / WindSp) * predicteds
      predicteds_old <- gm3toppm(predicteds_old, air_temp, air_press, 16.04)
      
      # Calculate log ratio
      
      log_prop_ratio_sigz <- (dgamma(sigz_prop, shape = 1.608356, rate = 0.7360556, log = TRUE) +
                                log_like(observeds, predicteds_prop, tau[i, ],
                                         data_size, classification, liketype, tau_WS, WindSp_hat)) -
                             (dgamma(sigz[i - 1], shape = 1.608356, rate = 0.7360556, log = TRUE) +
                                 log_like(observeds, predicteds_old, tau[i, ],
                                          data_size, classification, liketype, tau_WS, WindSp_hat))
      
      # Exp(log ratio)
      prop_ratio_sigz <- exp(log_prop_ratio_sigz)
      
      # Accept/reject
      if (prop_ratio_sigz >= 1){
        sigz[i] <- sigz_prop
        predicteds <- predicteds1
        acceptance_sigz[i] <- 1
      } else {
        uu <- runif(1)
        if (uu <= prop_ratio_sigz){
          sigz[i] <- sigz_prop
          predicteds <- predicteds1
          acceptance_sigz[i] <- 1
        } else {
          sigz[i] <- sigz[i - 1]
          predicteds <- predicteds
          acceptance_sigz[i] <- 0
        }
      }
    }
  
  
    ####### Q
    
    # Adapt
    if ((i %% 40 == 0 & i >= 50 & i < 1000) | (i == 50)){
        if (mean(acceptance[(i - 40):(i - 1)]) < 0.1){ # Adapt Q SD on proposal
        Qsd <- Qsd * 0.9
      } else if (mean(acceptance[(i - 40):(i - 1)]) > 0.8) {
        Qsd <- Qsd * 1.1
      }
    } else if (i %% 100 == 0 & i >= 1000 & i <= 20000){
      if (mean(acceptance[(i - 100):(i - 1)]) < 0.1){
        Qsd <- Qsd * 0.9
      } else if (mean(acceptance[(i - 100):(i - 1)]) > 0.8) {
        Qsd <- Qsd * 1.1
      }
      if (i == 20000){
        print(paste0("Q acceptance rate at completion of adaptation: ", 
                   mean(acceptance)), quote = FALSE)
      }
    }
    
    # Sample 
    
    prop <- rnorm(1, mean = Q[i - 1], sd = Qsd) # propose value
    
    # Update predicted vals with newest ACCEPTED Q and WS (for denominator of log ratio)
    predicteds2 <- (Q[i - 1] / WindSp) * predicteds
    predicteds2 <- gm3toppm(predicteds2, air_temp, air_press, 16.04) # convert from g/m^3 to ppm
    
    # Update predicted vals with newest Q and WS (for numerator of log ratio)
    predicteds_Q <- (prop / WindSp) * predicteds
    predicteds_Q <- gm3toppm(predicteds_Q, air_temp, air_press, 16.04) # convert from g/m^3 to ppm
    
    # Calculate log ratio
    log_prop_ratio <- (log_Q(prop, Q_prior_sd) + log_like(observeds, predicteds_Q, 
                        tau[i, ], data_size, classification,
                        liketype, tau_WS, WindSp_hat)) -
                      (log_Q(Q[i - 1], Q_prior_sd) + log_like(observeds, predicteds2, 
                        tau[i, ], data_size, classification,
                        liketype, tau_WS, WindSp_hat))
    
    # exp(log ratio)
    prop_ratio <- exp(log_prop_ratio)
    
    # Accept/reject
    if (prop_ratio >= 1){
      Q[i] <- prop
      acceptance[i] <- 1
    } else {
      uu <- runif(1)
      if (uu <= prop_ratio){
        Q[i] <- prop 
        acceptance[i] <- 1
      } else {
        Q[i] <- Q[i - 1]
        acceptance[i] <- 0
      }
    }
    
    if (i == iterates){ # print end time
      print(Sys.time(), quote = FALSE)
    }
    
  }

  ################################################################################

  # Get data together to return from function

  acc_rate_Q <- mean(acceptance[iterates - 100:iterates]) # Acceptance rate for last 100 iterates
  
  Q_burned <- Q[(burnin + 1):iterates] # remove burnin
  sigy_burned <- sigy[(burnin + 1):iterates] # remove burnin
  sigz_burned <- sigz[(burnin + 1):iterates] # remove burnin
  tau_burned <- as.matrix(tau[((burnin + 1):iterates), ]) # remove burnin
  
  thins <- seq(from = 1, to = iterates - burnin, by = thin) # create thinning vector
  Q_thinned <- Q_burned[-thins] # thin Q
  sigy_thinned <- sigy_burned[-thins] # thin omega_y
  sigz_thinned <- sigz_burned[-thins] # thin omega_z
  tau_thinned <- as.matrix(tau_burned[-thins, ]) # thin tau
  
  Q_sum <- c(mean(Q_thinned), median(Q_thinned), sd(Q_thinned), var(Q_thinned), quantile(Q_thinned, c(0.025, 0.975)), quantile(Q_thinned, c(0.01, 0.99))) # summary stats for Q
  sigy_sum <- c(mean(sigy_thinned), median(sigy_thinned), sd(sigy_thinned), var(sigy_thinned), quantile(sigy_thinned, c(0.025, 0.975)), quantile(sigy_thinned, c(0.025, 0.975))) # summary stats for omega
  sigz_sum <- c(mean(sigz_thinned), median(sigz_thinned), sd(sigz_thinned), var(sigz_thinned), quantile(sigz_thinned, c(0.025, 0.975)), quantile(sigz_thinned, c(0.01, 0.99))) # summary stats for omega

  
  results <- list(acc_rate_Q, Q_thinned, tau_thinned, Q_sum, sigz_thinned, sigz_sum, sigy_thinned, sigy_sum) # list to return
                  
  return(results)
}
