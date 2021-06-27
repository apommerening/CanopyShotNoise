# Prototype of a crown IBM.
# Ap - Goettingen, 07.10.2020. Modified on 28.10.2020. 
# 

rm(list = ls())									
options(digits = 6, width = 50)

# Packages
library(spatstat)
library(Rcpp)

codeFile <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/CanopyIBM/"
sourceCpp(paste(codeFile, "Movement.cpp", sep = ""))
sourceCpp(paste(codeFile, "GrowthInteraction.cpp", sep = ""))


# Settings 
set.seed(round(runif(1, min = 1, max = 10000)))
xmax <- 100 # Observation window defined by xmax and ymax
ymax <- 100

# Poisson process for initialisation, the points are the crown centres
myLambda <- 0.025 # Density
pattern <- rpoispp(lambda = myLambda, win = owin(c(0, xmax), c(0, ymax)), nsim = 1)
# Calculate local density for each point
xlambda <- density(pattern, at = "points") # sigma = bw.diggle, # bw.ppl

# Auxiliary graphs to check up on density-dependent marking model
par(mar = c(2, 4.5, 0.5, 0.5))
plot(xlambda,  2.8 - 40 * xlambda, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE)#, ylim = c(0, max(myMarks)), xlim = c(0, max(xlambda)))
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

par(mar = c(2, 4.5, 0.5, 0.5))
plot(xlambda,  6 - 55.317 * xlambda - 56.2 * max(xlambda), las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE)#, ylim = c(0, max(myMarks)), xlim = c(0, max(xlambda)))
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)

# Density-dependent marking, i.e. the larger local density is the smaller is the crown radius
myMarks <- c()
prandom <- runif(length(xlambda), min = 0, max = 1)
m <- median(xlambda)
for(i in 1 : length(xlambda)) {
  scale <- 4 - 91 * xlambda[i] # 2.8
  shape <- 7 - 50.317 * xlambda[i] - 1.2 * m # max(xlambda)
  Ts <- 1 + 1 / (100 *  xlambda[i])
  myMarks[i] <-  scale * ((Ts / scale)^shape - log(1 - prandom[i]))^(1 / shape) # scale[1]
}

# Auxiliary graph to check on density trend
par(mar = c(2, 3.3, 0.5, 0.5))
plot(xlambda, myMarks, las = 1, ylab = "", xlab = "", 
     cex = .9, col = "black", pch = 16, axes = FALSE, ylim = c(0, max(myMarks)), xlim = c(0, max(xlambda)))
axis(1, lwd = 2, cex.axis = 1.8) 
axis(2, las = 1, lwd = 2, cex.axis = 1.8)
box(lwd = 2)
cor(xlambda, myMarks)


# Visualisation of intial state of tree crown with density map in the background
pattern <- ppp(pattern$x, pattern$y, xrange = c(0, xmax), yrange = c(0, ymax), marks = myMarks)
par(mar = c(0, 0, 0, 0))									
# contour(density(pattern), 10, axes = FALSE, main = "", col = "black")
plot(density(pattern, sigma = 9, edge = TRUE), main = "") # sigma = bw.scott(pattern)[1]
plot(pattern, main = "", legend = F, bg = "grey", fg = "grey", markscale = 1.0, add = T) # bg = "red"

# Model parameters (currently just assumed)
a <- 1.10019 # Potential growth
b <- -0.037251 # Potential growth
abdn <- c(0.142, 1.135, 1.863, 2.394) # Growth interaction based on Haebel et al. (2019), DF
simYears <- 50  # Number of simulation years
mortThres <- 0.90 # Mortality threshold: Lowest relative growth rate (growth multiplier)

# i <- 1
for(i in 1 : simYears) {
  # Growth
  potGrowth <- a * myMarks^{b} # Potential growth
  redGrowth <- simulateRelativeGrowthRate(data.frame(x = pattern$x, y = pattern$y, mark = pattern$marks), 
                                        abdn, xmax, ymax, potGrowth)
  myMarks <- myMarks * redGrowth # Modifier
  # Mortality
  pattern <- pattern[redGrowth > mortThres, ] # Keep only alive ones
  myMarks <- myMarks[redGrowth > mortThres] # Keep only alive ones
  # Movement
  maxmov <- 0.3 # Maximum movement radius in m
  maxAttempts <- 100 # Number of attempts to find a better location within maxmov
  # Visualisation of new state
  newCoords <- moveCrowns(i, maxmov, maxAttempts, pattern$x, pattern$y, pattern$marks, xmax, ymax, 
             abdn, FALSE)
  pattern$x <- newCoords$x
  pattern$y <- newCoords$y
  pattern <- ppp(pattern$x, pattern$y, xrange = c(0, xmax), yrange = c(0, ymax), marks = myMarks)
  par(mar = c(0, 0, 0, 0))									
  plot(density(pattern, sigma = 9, edge = TRUE), main = "") # sigma = bw.scott(pattern)[1]
  plot(pattern, main = "", legend = F, bg = "grey", fg = "grey", markscale = 1.0, add = T) # bg = "red"
  # print(length(pattern$x))
}
