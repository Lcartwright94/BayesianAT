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


file_name <- "B1"
bg_type <- "5th-percentile"
exp_setup <- "April 23 -- June 7"
datafile <- "GA-inversion-data-1.csv" # where to get data from
no_of_insts <- 7
inst_numbers <- 1:7
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
