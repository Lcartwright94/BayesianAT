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
