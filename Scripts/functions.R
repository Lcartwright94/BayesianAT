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


###### Find Gamma params

Findalphabeta_gamma <- function(pars,p5,p95) {
  if(any(pars < 0)) {
    return(Inf)
  } else {
    return( sum((qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2]) - c(1/p95^2,1/p5^2))^2/c(1/p95^2,1/p5^2)))
  }
}


###### Set a, b, c, d parameters

stability_classes <- data.frame(
  c(0.17992784, 0.94470, 24.1670, 2.5334),
  c(0.14505676, 0.93198, 18.3330, 1.8096),
  c(0.11025100, 0.91465, 12.5000, 1.0857),
  c(0.08473887, 0.86974, 8.3330, 0.72382),
  c(0.07500510, 0.83660, 6.2500, 0.54287),
  c(0.05437017, 0.81558, 4.1667, 0.36191)
)

###### get stability class

Stab_class <- function(data){
  Linv <- 1 / data$L
  sapply(Linv, function(x){
    if (x <= -0.12){1}
    else if (-0.12 <= x & x < -0.06){2}
    else if (-0.06 <= x & x < -0.02){3}
    else if (-0.02 <= x & x < 0.02){4}
    else if (0.02 <= x & x < 0.07){5}
    else if (0.07 <= x){6}
  })
}


################ Gaussian plume without rotation (without Q / U) using matrix operations

gaussian_plume_mat <- function(X, Y, z, H, stability, sigy, sigz){
  concs <- matrix(nrow = nrow(X), ncol = ncol(X))
  sigma_y <- concs 
  sigma_z <- concs
  logX <- concs
  exp1 <- concs
  exp2 <- concs
  logX <- concs
  
  logX[X < 0] <- 0.01
  logX[X >= 0] <- log(X[X >= 0])
  
  concs[X < 0] <- 0
  
  if (dim(X)[2] == 1){ # points
      TH <- 0.01745 * (stability[, 3] - stability[, 4] * as.vector(logX) + stability[, 4] * log(1000))
      sigma_y <- sigy * as.vector(0.46511628 * X * tan(TH))
      #sigma_y <- as.vector(stability[, 3] * X ^ stability[, 4])
      sigma_z <- sigz * as.vector(stability[, 1] * X ^ stability[, 2])
      
      exp1 <- exp(-((as.vector(Y) ^ 2) / (2 * sigma_y ^ 2)))
      exp2 <- exp(-((z - H) ^ 2) / (2 * sigma_z ^ 2)) + exp(-((z + H) ^ 2) / (2 * sigma_z ^ 2))
      exp11 <- exp1 / (sigma_y * sqrt(2 * pi))
      exp22 <- exp2 / (sigma_z * sqrt(2 * pi))
      preds <- exp11 * exp22
      concs[X >= 0] <- preds[X >= 0]
      conc <- concs
  } else { # paths
    #ones <- matrix(rep(1, nrow(concs)), nrow = nrow(concs))
    TH <- 0.01745 * (stability[, 3] - stability[, 4] * logX + stability[, 4] * log(1000))
    #sigma_y <- sigy * 0.46511628 * X * tan(TH)
    sigma_y <- 0.46511628 * X * tan(TH)
    #sigma_y <- stability[, 3] * X ^ stability[, 4]
    sigma_z <- sigz * stability[, 1] * X ^ stability[, 2]
    #sigma_z <- stability[, 1] * X ^ stability[, 2]
    Z <- z %*% t(matrix(rep(1, ncol(X)), nrow = ncol(X)))
    
    exp1 <- exp(-((Y ^ 2) / (2 * sigma_y ^ 2)))
    exp2 <- exp(-((Z - H) ^ 2) / (2 * sigma_z ^ 2)) + exp(-((Z + H) ^ 2) / (2 * sigma_z ^ 2))
    exp11 <- exp1 / (sigma_y * sqrt(2 * pi))
    exp22 <- exp2 / (sigma_z * sqrt(2 * pi))
    preds <- exp11 * exp22
    concs[X >= 0] <- preds[X >= 0]
    conc <- rowMeans(concs)
  }
  return(as.vector(conc))
}


################ Convert to ppm

gm3toppm <- function(concs, temp, press, molar_mass){
  moles <- concs / molar_mass
  vol <- 1e6 * moles * 8.3144598 * temp / press
  return(vol)
}

########### Log likelihood (all obs, for Q)

log_like <- function(obs, preds, Tau, N, classification, liketype, tau_WS, WindSp_hat){
  if (length(levels(classification)) == 1){
    if (tau_WS == 1){
      TauVec <- Tau * (WindSp_hat ^ 4)
    } else {
      TauVec <- rep(as.numeric(Tau), N)
    }
  } else {
    if (tau_WS == 1){
      tmat <- model.matrix(~ classification + 0)
      TauVec <- as.vector(tmat %*% Tau)
      TauVec <- (WindSp_hat ^ 4) * TauVec
    } else {
      tmat <- model.matrix(~ classification + 0)
      TauVec <- as.vector(tmat %*% Tau)
    }
  }
  TauMatrix <- Diagonal(x = TauVec)

  if (liketype == "normal"){
    func <- sum((1 / 2) * log(TauVec)) - (N / 2) * log(2) - (N / 2) * log(pi) - 
      (1 / 2) * as.numeric(((t(obs - preds) %*% TauMatrix %*% (obs - preds))))
  } else if (liketype == "laplace"){
    func <- - N * log(2) + sum(log(TauVec)) - 
      sum(TauVec * abs(obs - preds))
  }
  return(func)
}

########### Log likelihood (single obs, in case needed)

log_like_sing <- function(obs, preds, Tau, classification, liketype, tau_WS, WindSp_hat){
  if (length(levels(classification)) == 1){
    if (tau_WS == 1){
      TauVec <- Tau * (WindSp_hat ^ 4)
    } else {
      TauVec <- rep(as.numeric(Tau), length(preds))
    }
  } else {
    if (tau_WS == 1){
      tmat <- model.matrix(~ classification + 0)
      TauVec <- as.vector(tmat %*% Tau)
      TauVec <- (WindSp_hat ^ 4) * TauVec
    } else {
      tmat <- model.matrix(~ classification + 0)
      TauVec <- as.vector(tmat %*% Tau)
    }
  }
  
  if (liketype == "normal"){
    func <- (1 / 2) * log(TauVec) - (1 / 2) * log(2) - (1 / 2) * log(pi) - 
      (TauVec / 2) * ((obs - preds) ^ 2)
  } else if (liketype == "laplace"){
    func <- log(TauVec) - log(2) - 
      (TauVec * abs(obs - preds))
  }
  return(func)
}


########## Log Likelihood for Tau

log_like_tau <- function(obs, preds, Tau, classification, sizes, tau_WS, WindSp_hat){
  if (length(levels(classification)) == 1){
    tmat <- matrix(rep(1, sizes), nrow = sizes, ncol = 1)
  } else {
    tmat <- model.matrix(~ classification + 0)
  }
  
  tmat_obs <- scale(t(tmat), center = FALSE, scale = 1 / obs)
  tmat_preds <- scale(t(tmat), center = FALSE, scale = (1 / (preds)))
  if (tau_WS == 1){
    tmat_WS <- t(scale(t(tmat), center = FALSE, scale = 1 / (WindSp_hat ^ 4)))
    tmat_WS_log <- tmat_WS
    tmat_WS_log[tmat_WS_log > 0] <- log(tmat_WS[tmat_WS > 0])
  }
  
  if (liketype == "normal"){
    if (tau_WS == 1){
      tmat_diffs <- (tmat_obs - tmat_preds) ^ 2
      tmat_diffs_WS <- t(scale(tmat_diffs, center = FALSE, scale = 1 / (WindSp_hat ^ 4)))
      diff_sums <- as.vector(colSums(tmat_diffs_WS))
    } else {
      tmat_diffs <- t((tmat_obs - tmat_preds) ^ 2)
      diff_sums <- as.vector(colSums(tmat_diffs))
    }
  } else if (liketype == "laplace"){
    if (tau_WS == 1){
      tmat_diffs <- abs(tmat_obs - tmat_preds)
      tmat_diffs_WS <- t(scale(tmat_diffs, center = FALSE, scale = 1 / (WindSp_hat ^ 4)))
      diff_sums <- as.vector(colSums(tmat_diffs_WS))
    } else {
      tmat_diffs <- t(abs(tmat_obs - tmat_preds))
      diff_sums <- as.vector(colSums(tmat_diffs))
    }
  }
  
  func <- c()
  func[Tau <= 0] <- 0
  
  if (liketype == "normal"){
    if (tau_WS == 1){
    func[Tau > 0] <- (sizes[Tau > 0] / 2) * log(Tau[Tau > 0]) - (sizes[Tau > 0] / 2) * log(2) - 
      (sizes[Tau > 0] / 2) * log(pi) + (1 / 2) * colSums(tmat_WS_log)[Tau > 0] -
      (Tau[Tau > 0] / 2) * diff_sums[Tau > 0]
    } else {
      func[Tau > 0] <- (sizes[Tau > 0] / 2) * log(Tau[Tau > 0]) - (sizes[Tau > 0] / 2) * log(2) - 
        (sizes[Tau > 0] / 2) * log(pi) - (Tau[Tau > 0] / 2) * diff_sums[Tau > 0]
    }
  } else if (liketype == "laplace"){
    if (tau_WS == 1){
      func[Tau > 0] <- -sizes[Tau > 0] * log(2) + sizes[Tau > 0] * log(Tau[Tau > 0]) +
        colSums(tmat_WS_log)[Tau > 0] - (Tau[Tau > 0] * diff_sums[Tau > 0])
    } else {
      func[Tau > 0] <- - sizes[Tau > 0] * log(2) + sizes[Tau > 0] * log(Tau[Tau > 0]) - 
        (Tau[Tau > 0] * diff_sums[Tau > 0])
    }
  }
  
  return(func)
}


########## Log Q prior

log_Q <- function(Q, scale){
  if (Q <= 0){
    func <- -Inf
  } else {
  func <- (1 / 2) * log(2 / pi) - log(scale) - 
    ((Q ^ 2) / (2 * (scale ^ 2)))
  }
  return(func)
}


########## Log tau prior

log_tau <- function(tau, alpha, beta, classification){
  func <- c()
  func[tau <= 0] <- -Inf
  
  func[tau > 0] <- alpha * log(beta) - log(gamma(alpha)) + 
    (alpha - 1) * log(tau[tau > 0]) - beta * tau[tau > 0]
  return(func)
}



