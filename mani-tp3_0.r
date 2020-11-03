## Description : Lab 3 of Manifold Learning 
## Date: November 2016 by jc

## 1. 3-sensor data set ####
## Artificial data as from P. Desmartines, PhD Tesis 1994
generateData <- function(n) {
  require(pdist)
  
  # these sensors where selected randomly
  sensors <- matrix(ncol = 3, data = 
    c(0.026, 0.236, -0.653, 0.310, 0.507, -0.270, -0.466,  -0.140, 0.353, -0.473,
      0.241, 0.193, 0.969, 0.094, 0.756, -0.978, -0.574, -0.502, -0.281, 0.993,
      0.026, -0.913, -0.700, 0.876, 0.216, -0.739, 0.556, -0.155, 0.431, 0.411))
  
  # draw random points on the 3d unit cube
  
  points <- matrix(ncol=3, data = runif(3*n,-1,1))
  
  # We describe each point as the distance to sensors : intrinsic dimension = 3
  # while extrinsic dimension = 10
  return (as.matrix(pdist(points, sensors)))
}

res100   <- generateData(100)
res1000  <- generateData(1000)
res10000 <- generateData(10000)

## 2. PCA estimator ####
pca_res100 <- prcomp(res100)
pca_res1000 <- prcomp(res1000) 
pca_res10000 <- prcomp(res10000)
summary(pca_res100)
## 3. Correlation Dimension Estimator ####
corrDim <- function(data, epsilon = 10^seq(-2, 1, length.out = 100)){
  n <- nrow(data)
  C <- numeric(length(epsilon))
  dn <- c()
  for (i in 1:(n-2)){
    dn <- c(dn, sqrt(colSums(apply(data[(i+1):n,], 1, function(j) (j-data[i,]))^2)))
  }
  dn <- c(dn, sqrt(sum((data[n]-data[n-1])^2)))
  for (e in 1:length(epsilon)){
    C[e] <- sum(dn <= epsilon[e])
  }
  C <- C*2/(n*(n-1))
  return (list(epsilon=epsilon, C=C))
}


derivate <- function(x, y) {
  (y[-1]-y[-length(y)])/(x[-1]-x[-length(x)])
}

X <- generateData(1000)
Xdim <- corrDim(X, epsilon = 10^seq(-4, 1, length.out = 100))
dev <- derivate(log10(Xdim$epsilon), log10(Xdim$C))
dev
# Plot C2 vs epsilon (in log-log)
plot(log10(Xdim$epsilon), log10(Xdim$C), type = 'l', 
     xlab = expression(log(epsilon)), ylab = expression(log(C(epsilon))))

# Plot d log(C2) vs d log(epsilon) 
dev <- derivate(log10(Xdim$epsilon), log10(Xdim$C))

plot(log10(Xdim$epsilon[which(dev>0)]), dev[which(dev>0)], 
     type = 'l')


Xdim100  <- corrDim(res100, epsilon = 10^seq(-2, 1, length.out = 100))
Xdim1000 <- corrDim(res1000, epsilon = 10^seq(-2, 1, length.out = 100))
#Xdim10000 &lt;- corrDim(res10000, epsilon = 10^seq(-4, 1, length.out = 100))

# Plot C2 vs epsilon (in log-log)
plot(log10(Xdim100$epsilon), log10(Xdim100$C), type = 'l', ylim = c(0, 6), 
     xlab = expression(log(epsilon)), ylab = expression(log(C(epsilon))))
lines(log10(Xdim1000$epsilon), log10(Xdim1000$C), col = 2,
     xlab = expression(log(epsilon)), ylab = expression(log(C(epsilon))))
 
# Plot d log(C2) vs d log(epsilon) 
plot(log10(Xdim100$epsilon[-1]), 
     derivate(log10(Xdim100$epsilon), log10(Xdim100$C)), 
     type = 'l')
lines(log10(Xdim1000$epsilon[-1]), 
      derivate(log10(Xdim1000$epsilon), log10(Xdim1000$C)), 
      col = 2)

# 4. Spiral
...

</pre></body></html>Ztext/plainUUTF-8_Vhttps://moodle.univ-lyon2.fr/pluginfile.php/348034/mod_resource/content/0/mani-tp3_0.rP