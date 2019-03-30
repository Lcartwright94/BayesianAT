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


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(fdrtool))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(parallel))

##############################################################

### Read in functions needed for background estimation

source("upwind-downwind.R")
source("5th-percentile-background.R")
source("functions.R")
source("MCMC.R")

##############################################################

### Values you may wish to change

sigy_track <- 1 # Learn about sigma_y scalar
sigz_track <- 1 # Learn about sigma_z scalar
stab_track <- 1 # Track uncertainty by stability class in MCMC
Tau_track <- 1 # track tau, also track tau by inst type
tau_WS <- 1 # Create WS_hat variable

liketype <- "normal" # Form of likelihood function (laplace or normal)

wind_speed_cutoff <- 0 # Remove all observations with wind speed below this cutoff before running MCMC
samples <- 100 # how many pieces to cut line integral into
angle <- 45 # Upwind cutoff angle

iterates <- 60000 # How many MCMC samples
burnin <- 20000 # How many to throw away as burnin
thin <- 10 # How often to thin

##############################################################

params <- c("params-B1.R", "params-F1.R", "params-E1.R", "params-P1.R", "params-BFEP1.R", "params-B2.R", "params-F2.R", "params-E2.R", "params-P2.R", "params-BFEP2.R")
cores <- 4

##############################################################



#---------------------------------------------------------------------------------------------------#
#####################################################################################################
#####################################################################################################
########################################## START AUTO ###############################################
#####################################################################################################
#####################################################################################################




mclapply(1:length(params), function(p){

# Read in parameters

source(params[p])

print(paste0("Starting ",file_name), quote = FALSE)
  
################################################################

### Read in data file and begin background estimation

inv_data <- read.csv(datafile, header = TRUE)
inv_data <- filter(inv_data, wind_speed >= wind_speed_cutoff)

MyDate <- as.POSIXct(inv_data$Time, format = "%Y-%m-%d %H:%M")
Month <- month(MyDate)
Day <- day(MyDate)
Year <- year(MyDate)
Hour <- hour(MyDate)
Min <- minute(MyDate)
MyTime <- as.POSIXct(paste0(Year,"/",Month,"/",Day," ",Hour,":",Min))
inv_data$Time <- MyTime

# Filter out times when source is inactive
inv_data <- filter(inv_data, release_rate >= source_on_rate - 0.1 & release_rate < source_on_rate + 0.2)

# Filter out unwanted instruments
if (length(inst_numbers) == 1){
  inv_data <- filter(inv_data, inst_no == inst_numbers)
} else {
  inst_remove <- c()
  for (i in 1:nrow(inv_data)){
    if (length(which(inst_numbers == inv_data$inst_no[i])) > 0){
      inst_remove[i] <- 0
    } else {
      inst_remove[i] <- i
    }
  }
  inst_remove <- inst_remove[inst_remove != 0]
  if (length(inst_remove) > 0){
    inv_data <- inv_data[-inst_remove, ]
  }
}

# Stop running if all observations have been removed, warn if not many observations left
if (nrow(inv_data) == 0){
    stop("No data values!")
    terminate
} else if (nrow(inv_data) <= 50){
    print(paste0("WARNING: low number of observations (", nrow(inv_data), ")"), quote = FALSE)
}

### Stability Class
inv_data$Stability <- Stab_class(inv_data)

### Add plume cutoff angle

inv_data$angle <- rep(angle, nrow(inv_data))

### Background

inv_data_temp <- inv_data
inv_data <- Upwind.Downwind(inv_data = inv_data_temp) # sort into upwind/downwind observations

upwind_count <- nrow(inv_data[inv_data$Position == "upwind", ])
downwind_count <- nrow(inv_data[inv_data$Position == "downwind", ])
# Add warning if upwind count low, stop running if upwind count is 0
if (upwind_count < 1){
    stop("No upwind values!")
    terminate
} else if (upwind_count <= 40){
    print(paste0("WARNING: Upwind count low (", upwind_count, ")"), quote = FALSE)
}

# Run the background estimation function
inv_data_temp <- inv_data
inv_data <- fifth.perc.bg(inv_data = inv_data_temp)


######################################################################

# Add upwind and downwind counts

upwind_total <- rep(upwind_count, nrow(inv_data))
downwind_total <- rep(downwind_count, nrow(inv_data))

# Add inst group

inst_group <- inv_data$inst_group

######################################################################

# create data frame for the MCMC

mcmc_data <- data.frame(inv_data[, c("air_temp", "air_pressure",
                                 "wind_speed", "wind_dir", "L", "inst_no",
                                 "source_x", "source_y", "z", "rot_x",
                                 "rot_y", "rot_las_x", "rot_las_y", "bg_conc", "bg_sd", "smalln",
                                 "release_rate", "Stability", "inst_name")], upwind_total, downwind_total,
                                  inst_group)
names(mcmc_data) <- c("air_temp", "air_pressure",
                  "wind_speed", "wind_dir", "L", "inst_no",
                  "source_x", "source_y", "z", "x1",
                  "y1", "x2", "y2", "bg_conc", "bg_sd", "smalln",
                  "release_rate", "Stability", "inst_name", "upwind_total", "downwind_total",
                  "inst_group")

################################################################################

#### Run MCMC

results <- mcmc_run(mcmc_data, H, samples, iterates, burnin, thin)

# ORDER in "results" list:

#acc_rate_Q, Q_thinned, tau_thinned, Q_sum, omegaz_thinned, omegaz_sum, omegay_thinned, omegay_sum



################################################################################

traces <- data.frame(results[[2]], results[[7]], results[[5]]) # Q, omega_y, omega_z
names(traces) <- c("Q", "omega_y", "omega_z")
traces2 <- as.mcmc.list(lapply(traces, mcmc))

Q <- traces2[[1]]
omega_y <- traces2[[2]]
omega_z <- traces2[[3]]

if (exp_setup == "April 23 -- June 7"){
  true_Q <- 0.09667
} else if (exp_setup == "June 8 -- June 12"){
  true_Q <- 0.08333
} else if (true_val_known == "No"){
  true_Q <- mean(Q)
}

################################################################################

# Generate plots

png(file = paste0(file_name,"-plots.png"))
par(mfrow = c(3, 3))

# Q
hist(Q, main = paste0("Histogram of Q, N = ", nrow(mcmc_data)), freq = FALSE, xlab = "Q", col = "skyblue", xlim = c(min(traces$Q, true_Q), max(traces$Q, true_Q)))
if (exp_setup == "April 23 -- June 7" || exp_setup == "June 8 -- June 12" || (exp_setup == "Other" & true_val_known == "Yes")){
  abline(v = true_Q, col = "red", lty = "dashed", lwd = 2)
}
plot(Q, density = FALSE, auto.layout = FALSE, col = "indianred1", main = "Trace of Q", ylab = "Q")
autocorr.plot(Q, auto.layout = FALSE, main = "Autocorrelation of Q")

# omega_y
hist(omega_y, main = "Histogram of omega_y", freq = FALSE, col = "skyblue", xlab = "omega_y")
plot(omega_y, density = FALSE, auto.layout = FALSE, col = "indianred1", main = "Trace of omega_y", ylab = "omega_y")
autocorr.plot(omega_y, auto.layout = FALSE, main = "Autocorrelation of omega_y")

# omega_z
hist(omega_z, main = "Histogram of omega_z", freq = FALSE, col = "skyblue", xlab = "omega_z")
plot(omega_z, density = FALSE, auto.layout = FALSE, col = "indianred1", main = "Trace of omega_z", ylab = "omega_z")
autocorr.plot(omega_z, auto.layout = FALSE, main = "Autocorrelation of omega_z")
dev.off()


### Create .csv file of summary stats

summaries <- rbind(c("parameter", "mean", "median", "sd", "var", "lower95", "upper95", "lower99", "upper99"),
                        c("Q", results[[4]]),
                        c("omega_y", results[[8]]),
                        c("omega_z", results[[6]]),
                        c(" ", " ", " ", " ", " ", " ", " ", " ", " "),
                        c("Upwind Count: ", upwind_count, "Downwind Count:", downwind_count, " ", " ", " ", " ", " "),
                        c("Q acceptance rate:", results[[1]], " ", " ", " ", " ", " ", " ", " "))


# write csv
write.csv(summaries, file = paste0(file_name,"-summary.csv"), row.names = FALSE)

# save entire list
save(results, file = paste0(file_name,"-res.RData"))

print(paste0("Completed ",file_name), quote = FALSE)

}, mc.cores = cores)

print("All finished!", quote = FALSE)
