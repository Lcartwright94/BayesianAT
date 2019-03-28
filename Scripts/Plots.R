library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)

source("upwind-downwind.R")
source("functions.R")
source("5th-percentile-background.R")

H <- 0.3
insts <- c("EC.A", "EC.C", "EC.D", "EC.E", "Picarro.East", "Picarro.West", "P1", "P2", "P3", "P4", "P5", "P6", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13")
names <- c("B1", "F1", "E1", "P1", "BFEP1", "B2", "F2", "E2", "P2", "BFEP2")

red <- "#FD7F7F"
blue <- "#7B84FC"

par(mfrow = c(1, 1))

### Script containing code to reproduce plots in AMT paper

# Setup plots

source_x <- -21.78
source_y <- 21.09
# Add values of x1, y1, x2 and y2 (location of lasers, reflectors, towers... whatever applicable)
x1 <- c(-27.44, -10.29, 3.48, 11.10, 6.96, 2.04, -2.93, 
        -10.28, -29.82, -5.80, 20.59, -19.38, 6.98,
        -53.38, -1.88, 35.12, -1.88, 35.12, 18.62,
        0, -30.62, -8.91, 7.31, 
        31.15, -71.69) #reflector/tower location
y1 <- c(52.54, 41.34, 30.66, 18.16, 2.72, -12.17, -22.27, 
        41.30, 19.70, 23.00, 29.69, 10.54, 2.72,
        62.13, 61.63, 42.63, 61.63, 42.63, -25.87,
        0, 66.73, -19.88, -12.93, 
        -64.81, 53.84) #reflector/tower location
x2 <- c(rep(-45.75, 13), rep(-29.12, 3), rep(-67.12, 3), rep(NA, 6)) # laser location (NA for towers) #-45.75
y2 <- c(rep(16.20, 13), rep(-39.71, 3), rep(27.22, 3), rep(NA, 6)) # laser location (NA for towers) #16.20

par(mfrow = c(2, 1))

# 5.8 g/min setup

plot(source_x, source_y, pch = 16, col = "black", xlim = c(-70, 40), ylim = c(-65, 65), 
     xlab = "x (m)", ylab = "y (m)", cex.axis = 1.2, cex.lab = 1.2) #source

points(x1[20:23], y1[20:23], pch = 8, col = "blue") #EC towers
points(x1[24:25], y1[24:25], pch = 8, col = "aquamarine4") #Picarro analysers
points(x2[1], y2[1], pch = 17, col = "purple") #boreal laser
points(x2[c(14, 17)], y2[c(14, 17)], pch = 17, col = "red") #FTIR spectrometers
points(x1[14:19], y1[14:19], pch = 18, col = "red") #FTIR reflectors

points(x1[1:7], y1[1:7], pch = 18, col = "purple") #boreal may reflectors

lines(c(x2[1], x1[1]), c(y2[1], y1[1]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[2]), c(y2[1], y1[2]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[3]), c(y2[1], y1[3]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[4]), c(y2[1], y1[4]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[5]), c(y2[1], y1[5]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[6]), c(y2[1], y1[6]), lty = "dashed", col = "purple") #boreal may
lines(c(x2[1], x1[7]), c(y2[1], y1[7]), lty = "dashed", col = "purple") #boreal may

lines(c(x2[14], x1[14]), c(y2[14], y1[14]), lty = "dashed", col = "red") #FTIR
lines(c(x2[14], x1[15]), c(y2[14], y1[15]), lty = "dashed", col = "red") #FTIR
lines(c(x2[14], x1[16]), c(y2[14], y1[16]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[17]), c(y2[17], y1[17]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[18]), c(y2[17], y1[18]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[19]), c(y2[17], y1[19]), lty = "dashed", col = "red") #FTIR

legend(-75, -10, legend=c("Boreal laser", "Boreal reflector", "FTIR spectrometer", "FTIR retro-reflector",
                          "EC tower", "Picarro analyser", "Source"),
       col = c("purple", "purple", "red", "red", "blue", "aquamarine4", "black"), 
       pch = c(17, 18, 17, 18, 8, 8, 16), cex = 1, bg = "transparent", bty = "n", x.intersp = 0.4, y.intersp = 0.6)


#5.0 g/min setup

plot(source_x, source_y, pch = 16, col = "black", xlim = c(-70, 40), ylim = c(-65, 65), 
     xlab = "x (m)", ylab = "y (m)", cex.axis = 1.2, cex.lab = 1.2) #source

points(x1[20:23], y1[20:23], pch = 8, col = "blue") #EC towers
points(x1[24:25], y1[24:25], pch = 8, col = "aquamarine4") #Picarro analysers
points(x2[1], y2[1], pch = 17, col = "purple") #boreal laser
points(x2[c(14, 17)], y2[c(14, 17)], pch = 17, col = "red") #FTIR spectrometers
points(x1[14:19], y1[14:19], pch = 18, col = "red") #FTIR reflectors

points(x1[8:13], y1[8:13], pch = 18, col = "purple") #boreal june reflectors

lines(c(x2[1], x1[8]), c(y2[1], y1[8]), lty = "dashed", col = "purple") #boreal june
lines(c(x2[1], x1[9]), c(y2[1], y1[9]), lty = "dashed", col = "purple") #boreal june
lines(c(x2[1], x1[10]), c(y2[1], y1[10]), lty = "dashed", col = "purple") #boreal june
lines(c(x2[1], x1[11]), c(y2[1], y1[11]), lty = "dashed", col = "purple") #boreal june
lines(c(x2[1], x1[12]), c(y2[1], y1[12]), lty = "dashed", col = "purple") #boreal june
lines(c(x2[1], x1[13]), c(y2[1], y1[13]), lty = "dashed", col = "purple") #boreal june

lines(c(x2[14], x1[14]), c(y2[14], y1[14]), lty = "dashed", col = "red") #FTIR
lines(c(x2[14], x1[15]), c(y2[14], y1[15]), lty = "dashed", col = "red") #FTIR
lines(c(x2[14], x1[16]), c(y2[14], y1[16]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[17]), c(y2[17], y1[17]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[18]), c(y2[17], y1[18]), lty = "dashed", col = "red") #FTIR
lines(c(x2[17], x1[19]), c(y2[17], y1[19]), lty = "dashed", col = "red") #FTIR

legend(-75, -10, legend=c("Boreal laser", "Boreal reflector", "FTIR spectrometer", "FTIR retro-reflector",
                          "EC tower", "Picarro analyser", "Source"),
       col = c("purple", "purple", "red", "red", "blue", "aquamarine4", "black"), 
       pch = c(17, 18, 17, 18, 8, 8, 16), cex = 1, bg = "transparent", bty = "n", x.intersp = 0.4, y.intersp = 0.6)

# Labels were added manually


##################

# Scaling sigma_y plots

mcmc_data <- read.csv("GA-inversion-data-1.csv", header = TRUE)
mcmc_data <- filter(mcmc_data, release_rate > 5.7 & release_rate < 6)
gam <- c(1, 2.5, 4)

mcmc_data$angle <- rep(45, nrow(mcmc_data))
mcmc_data$Stability <- Stab_class(mcmc_data)

inst <- 20

mcmc_data <- filter(mcmc_data, inst_no == inst)

mcmc_data_temp <- mcmc_data
mcmc_data <- Upwind.Downwind(mcmc_data_temp)

mcmc_data_temp <- mcmc_data
mcmc_data <- fifth.perc.bg(mcmc_data_temp)

data_size <- nrow(mcmc_data)
stability <- matrix(nrow = data_size, ncol = 4)
for (k in 1:data_size){
  stability[k, ] <- stability_classes[ , mcmc_data$Stability[k]]
}

# Calculate predicted values
predicteds <- matrix(nrow = data_size, ncol = 3)

par(mfrow = c(1, 3), mar=c(5, 5, 3, 3))

for (j in 1:3){
  preds <- gaussian_plume_mat(as.matrix(mcmc_data$rot_x), as.matrix(mcmc_data$rot_y), mcmc_data$z, H, stability, gam[j], 1)
  predicteds[, j] <- preds
  
  predicteds[, j] <- (0.09667 / mcmc_data$wind_speed) * predicteds[, j]
  predicteds[, j] <- gm3toppm(predicteds[, j], mcmc_data$air_temp, mcmc_data$air_pressure, 16.04)
  
  
  SE <- (predicteds[, j] - mcmc_data$bg_conc) ^ 2
  MSE <- round(sum(SE) / nrow(mcmc_data), 4)
  
  plot(mcmc_data$wind_dir, mcmc_data$bg_conc, xlim = c(280, 350), ylim = c(0, 2), col = alpha("red", 0.4),
       main = paste0("Scale = ",gam[j],"  MSE = ",MSE), xlab = "Wind direction (degrees)", ylab = "Concentration (ppm)", 
       cex.axis = 1.5, cex.lab = 1.5)
  points(mcmc_data$wind_dir, predicteds[, j], xlim = c(280, 350), ylim = c(0, 2), col = alpha("blue", 0.4))
}

##################

# Background subtracted plots

inv_data1 <- read.csv("GA-inversion-data-1.csv", header = TRUE)
inv_data2 <- read.csv("GA-inversion-data-2.csv", header = TRUE)
inv_data <- rbind(inv_data1, inv_data2)
inv_data$Time <- as.POSIXct(inv_data$Time, format = "%Y-%m-%d %H:%M:%S")
inv_data$Date <- inv_data$Time
inv_data$Instrument <- inv_data$inst_name
inv_data$Instrument <- as.character(inv_data$Instrument)
inv_data$Instrument <- factor(inv_data$Instrument, levels = insts)

inv_data$angle <- rep(45, nrow(inv_data))

inv_data_temp <- inv_data
inv_data <- Upwind.Downwind(inv_data_temp)

inv_data_temp <- inv_data
inv_data <- fifth.perc.bg(inv_data_temp)

p1 <- ggplot(inv_data) + geom_point(aes(x = wind_dir, y = Concentration, 
                                        colour = Instrument), alpha = 0.5) + theme_bw() + 
  xlab("Wind direction (degrees)") + ylab("Concentration (ppm)") +
  theme(axis.text.x= element_text(size = 12), axis.text.y= element_text(size = 12),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12)) + 
  ylim(-0.25, 13) + theme(plot.margin = unit(c(1,1,1,1), "lines"))


p2 <- ggplot(inv_data) + geom_point(aes(x = wind_dir, y = bg_conc, 
                                        colour = Instrument), alpha = 0.5) + theme_bw() + 
  xlab("Wind direction (degrees)") + ylab("Enhancement (ppm)") +
  theme(axis.text.x= element_text(size = 12), axis.text.y= element_text(size = 12),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  ylim(-0.25, 13) + theme(plot.margin = unit(c(1,1,1,2), "lines"))


ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")



##################

# Results plot (intervals)
# This will only work after results have been generated
# Will not produce "source off" plots

load(paste0(names[1],"-res.RData"))
B1 <- results[[8]]
load(paste0(names[2],"-res.RData"))
F1 <- results[[8]]
load(paste0(names[3],"-res.RData"))
E1 <- results[[8]]
load(paste0(names[4],"-res.RData"))
P1 <- results[[8]]
load(paste0(names[5],"-res.RData"))
BFEP1 <- results[[8]]
load(paste0(names[6],"-res.RData"))
B2 <- results[[8]]
load(paste0(names[7],"-res.RData"))
F2 <- results[[8]]
load(paste0(names[8],"-res.RData"))
E2 <- results[[8]]
load(paste0(names[9],"-res.RData"))
P2 <- results[[8]]
load(paste0(names[10],"-res.RData"))
BFEP2 <- results[[8]]


Qs <- data.frame(B1, F1, E1, P1, BFEP1, B2, F2, E2, P2, BFEP2)
Qs <- gather(Qs, key = "Instruments", "Q")
Qs$Q <- Qs$Q * 60
Qs$Period <- factor(c(rep("5.8", 360000/2), rep("5.0", 360000/2)), levels = c("5.8", "5.0"))
Qs$Instruments <- factor(Qs$Instruments, levels = c("BFEP2", "P2", "E2", "F2", "B2", "BFEP1", "P1", "E1", "F1", "B1"))

q1 <- function(x) {
  y = median(x)
  ymin <- quantile(x, 0.025)
  ymax <- quantile(x, 0.975)
  return(c(y = y, ymin = ymin[[1]], ymax = ymax[[1]]))
}

ggplot(Qs, aes(x = Instruments, y = Q, fill = Period, border = 0.3)) + geom_hline(yintercept = 5.8, col = red, linetype = "dashed", size = 0.8) + 
  geom_hline(yintercept = 5, col = blue, linetype = "dashed", size = 0.8) + geom_violin() + theme_bw() +
  coord_flip() + ylim(0, 8) + stat_summary(fun.y=median, geom="point", size=1, color="black") +
  stat_summary(fun.data=q1, geom="errorbar", width = 0.5, size = 0.3) + scale_fill_manual(values = c(red, blue)) +
  geom_vline(xintercept = c(1:10), size = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       axis.text.x= element_text(size = 12), axis.text.y= element_text(size = 12),
                                                       axis.title.x= element_text(size = 12), axis.title.y= element_text(size = 12),
                                                       legend.direction = "horizontal", legend.position = "top",
                                                       legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  labs(fill = expression(paste("Release period (g min"^"-1",")"))) + ylab(expression(paste("Q (g min"^"-1",")"))) + xlab(" ")



