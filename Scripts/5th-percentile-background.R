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


fifth.perc.bg <- function(inv_data){
  
  print("Starting background estimation", quote = FALSE)
  
  for (i in 1:nrow(inv_data)){
      temp <- filter(inv_data, inst_no == inv_data$inst_no[i])
      quant <- quantile(temp$Concentration, 0.05)
      inv_data$bg_conc[i] <- inv_data$Concentration[i] - quant[[1]]
      inv_data$smalln[i] <- 1
      inv_data$bg_sd[i] <- 1
  }
  
  print("Background estimation complete", quote = FALSE)
  
  return(inv_data)
  
}
