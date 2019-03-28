file_name <- "F2"
bg_type <- "5th-percentile"
exp_setup <- "June 8 -- June 12"
datafile <- "GA-inversion-data-2.csv" # where to get data from
no_of_insts <- 6
inst_numbers <- 14:19
tot_insts <- 25 # total number of instruments possible to use
true_val_known <- "Yes"
true_Q <- 0.08333
H <- 0.3
molar_mass <- 16.04

if (true_val_known == "No"){
  source_on_rate <- 0
} else {
  source_on_rate <- (true_Q * 60)
} 
