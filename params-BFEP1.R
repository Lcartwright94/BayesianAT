file_name <- "BFEP1"
bg_type <- "5th-percentile"
exp_setup <- "April 23 -- June 7"
datafile <- "GA-inversion-data-1.csv" # where to get data from
no_of_insts <- 19
inst_numbers <- c(1:7, 14:25)
tot_insts <- 25 # total number of instruments possible to use
true_val_known <- "Yes"
true_Q <- 0.09667
H <- 0.3
molar_mass <- 16.04

if (true_val_known == "No"){
  source_on_rate <- 0
} else {
  source_on_rate <- (true_Q * 60)
} 
