#### Function to sort observations into upwind or downwind

Upwind.Downwind <- function(inv_data){
  dim2 <- as.numeric(nrow(inv_data))

  for (i in 1:dim2){
    # If source is off, upwind immediately
    if (inv_data$release_rate[i] == 0){
      inv_data$Position[i] <- "upwind"
      inv_data$rot_x[i] <- NA
      inv_data$rot_y[i] <- NA
      inv_data$rot_las_x[i] <- NA
      inv_data$rot_las_y[i] <- NA
      inv_data$rot_angle[i] <- NA
      inv_data$rot_las_angle[i] <- NA
    } else {
      mm <- inv_data$inst_no[i] # set instrument number
      xlassub <- inv_data$x2[i] - inv_data$source_x[i]
      ylassub <- inv_data$y2[i] - inv_data$source_y[i] 
      xlocsub <- inv_data$x1[i] - inv_data$source_x[i] 
      ylocsub <- inv_data$y1[i] - inv_data$source_y[i] # subtract off source
      theta <- ((inv_data$wind_dir[i] - 270) %% 360) * pi / 180 # convert to radians
      xyrota <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2, byrow = TRUE)
      if (is.na(inv_data$x2[i]) == TRUE){ # rotated coordinates for lasers
        inv_data$rot_las_x[i] <- NA 
        inv_data$rot_las_y[i] <- NA
      } else { 
        rot_laser <- round(xyrota %*% c(xlassub, ylassub), 5)
        inv_data$rot_las_x[i] <- rot_laser[1]
        inv_data$rot_las_y[i] <- rot_laser[2]
      } #rotated coordinates for reflectors/towers
      rotvec <- round(xyrota %*% c(xlocsub, ylocsub), 5) # this way matches python's .dot(vec) thing
      inv_data$rot_x[i] <- rotvec[1]
      inv_data$rot_y[i] <- rotvec[2]
      # Find rotated angles & assign wind position
      # rotate reflectors/tower
      rot_x <- inv_data$rot_x[i]
      rot_y <- inv_data$rot_y[i]
      if (rot_x < 0 & rot_y < 0){
        inv_data$rot_angle[i] <- atan(rot_y / rot_x) - round(pi, 5)
      } else if (rot_x < 0 & rot_y >= 0){
        inv_data$rot_angle[i] <- atan(rot_y / rot_x) + round(pi, 5)
      } else {
        inv_data$rot_angle[i] <- atan(rot_y / rot_x)
      }
      inv_data$rot_angle[i] <- round(round(inv_data$rot_angle[i], 5) * 180 / pi, 5) # convert to degrees
      # now rotate laser
      if (is.na(inv_data$x2[i]) == TRUE){
        inv_data$rot_las_angle[i] <- NA
      } else {
        rot_x <- inv_data$rot_las_x[i]
        rot_y <- inv_data$rot_las_y[i]
        if (rot_x < 0 & rot_y < 0){
            inv_data$rot_las_angle[i] <- atan(rot_y / rot_x) - round(pi, 5) # inverse tan automatically returns values in quadrants 1 and 4, so need to adjust for quadrants 2 and 3
        } else if (rot_x < 0 & rot_y >= 0){
          inv_data$rot_las_angle[i] <- atan(rot_y / rot_x) + round(pi, 5)
        } else {
          inv_data$rot_las_angle[i] <- atan(rot_y / rot_x)
        }
        inv_data$rot_las_angle[i] <- round(round(inv_data$rot_las_angle[i], 5) * 180 / pi, 5) # convert to degrees
      } 
      
      # assign wind position
      if (is.na(inv_data$x2[i]) == TRUE){ # Towers
        if (inv_data$rot_x[i] < 0 || inv_data$rot_angle[i] < -inv_data$angle[i] || inv_data$rot_angle[i] > inv_data$angle[i]){
          inv_data$Position[i] <- "upwind"
        } else {
          inv_data$Position[i] <- "downwind"
        } # if both end points have same sign on angle, just check if an end point is within +-angle degrees
      } else if ((inv_data$rot_angle[i] < 0 & inv_data$rot_las_angle[i] < 0) || (inv_data$rot_angle[i] >= 0 & inv_data$rot_las_angle[i] >= 0)) {
        if ((inv_data$rot_angle[i] >= -inv_data$angle[i] & inv_data$rot_angle[i] <= inv_data$angle[i]) || (inv_data$rot_las_angle[i] >= -inv_data$angle[i] & inv_data$rot_las_angle[i] <= inv_data$angle[i])){
          inv_data$Position[i] <- "downwind"
        } else {
          inv_data$Position[i] <- "upwind"
        } # If either end point within +-angle degrees
      } else if ((inv_data$rot_angle[i] >= -inv_data$angle[i] & inv_data$rot_angle[i] <= inv_data$angle[i]) || (inv_data$rot_las_angle[i] >= -inv_data$angle[i] & inv_data$rot_las_angle[i] <= inv_data$angle[i])){ 
        inv_data$Position[i] <- "downwind"
      } else { # If both end points outside +-angle degrees, check x-intercept
        x <- c(inv_data$rot_x[i], inv_data$rot_las_x[i])
        y <- c(inv_data$rot_y[i], inv_data$rot_las_y[i])
        xint <- x[1] - ((x[2] - x[1]) / (y[2] - y[1]) * y[1]) # x-intercept of a straight line is x_1 - y_1 / m, where m is gradient
        if (xint < 0){
          inv_data$Position[i] <- "upwind"
        } else {
          inv_data$Position[i] <- "downwind"
        }
      }
    }
  }
  return(inv_data)
}
