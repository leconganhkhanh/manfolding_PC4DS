## File: mani-tp1.r
## Description: Solution du TP1 cours modèles classificatoires
## Date : jan 2017

# [Ex. 1] Simulation par la méthode de rejet ######################

# 1. Fonction de densité d'une va triagulaire
# Input
#    x :  abscisse
f_tri <- function(x) {
  SUP01 <- (0 <= x) & (x < 1)
  SUP12 <- (1 <= x) & (x < 2) 
  if (SUP01) {
    return(x)
  } else {if (SUP12) { 
    return(2 - x)
  } else {
    return(0)
  }
  }
}

# Version optimizée
# f_tri2 <- function(x) {
#   SUP01 <- (0 <  x) & (x < 1)
#   SUP12 <- (1 <= x) & (x < 2) 
#   ifelse(SUP01, x, ifelse(SUP12, 2 - x, 0))
# }

# Test f_tri
# x <- seq(-1, 3, length.out = 101)
# f_tri_x <- numeric(length(x)) 
# for(i in seq_along(x)) f_tri_x[i] <- f_tri(x[i])

#plot(x, f_tri_x, type = 'l')
#rug(Observations)


# 2. 
# Input
#    fx   :  fonction de densité
#    a, b :  bornes du support de fx
#    M    :  borne supérieure de fx
rejection <- function(fx, a, b, M) {
  while (TRUE) { 
    x <- runif(1, a, b)
    y <- runif(1, 0, M)  
    if (y < fx(x)) return(x)
  }
}


# 3. Test rejection
nreps <- 1000
Observations <- numeric(nreps)
for (i in seq_along(Observations)) 
  Observations[i] <- rejection(f_tri, 0, 2, 1)


# [Ex. 2] Estimation de la densité par histogramme ##################

# 1. 
nreps <- 1000
Observations <- numeric(nreps)
for (i in seq_along(Observations)) 
  Observations[i] <- rejection(f_tri, 0, 2, 1)

# 2. 
m1 <- 1 + 10 * log(nreps) / 3
m2 <- sqrt(nreps)

layout(matrix(1:5, 1, 5))
for (b in c(5, 25, 500, m1, m2))
  hist(Observations, breaks = seq(0, 2, length.out = b + 1))

# 3. Calibrate hist (assumes that density's support is unit interval)
riskhist <- function(obs, m, xlim = c(0, 1)) {
  obs01  <- (obs - xlim[1]) / (xlim[2] - xlim[1])
  h      <- 1 / m
  n      <- length(obs)
  breaks <- seq(0, 1, length.out = m + 1)
  p_hat  <- hist(obs01, plot = FALSE, breaks = breaks)$counts / n
  
  res <- 2 / h / (n - 1) - (n + 1) / (n - 1) / h * sum(p_hat^2)
  return(res) #(m * sum(p_hat^2))
}

Mgrid <- 21:100
J     <- numeric(length(Mgrid))
for (m in seq_along(Mgrid)) J[m] <- riskhist(Observations, Mgrid[m], xlim = c(0, 2))
layout(1)
plot(1 / Mgrid, J, type = 'l')


# [Ex. 3] Estimation de la densité par noyau ##################

# 1.
nreps <- 1000
Observations <- numeric(nreps)
for (i in seq_along(Observations))
  Observations[i] <- rejection(f_tri, 0, 2, 1)

# 2.
s <- min(sd(Observations), diff(quantile(Observations, c(1/4, 3/4))) / 1.34)
h <- 1.06 * s / nreps ^ (1/5)

layout(matrix(1:4, ncol = 4))
for (b in c(0.01, 0.1, 1, h)) {
  plot(density(Observations, bw = b, from = 0, to = 2),
       main = paste("fenêtre =", b))
  rug(Observations)
}

# test the different kernels
#layout(matrix(1:8, nrow = 2))
#kernels <- c("gaussian", "epanechnikov", "rectangular",
#             "triangular", "biweight", "cosine", "optcosine")
#for(k in kernels){
#  plot(density(Observations, kernel = k, from = 0, to = 2),
#       main = paste("noyau =", k))
#  rug(Observations)
#}


# 3. Calibrate KDE
Kast <- function(x) dnorm(x, sd = 2) - 2 * dnorm(x)

riskKDE <- function(obs, h, xlim = c(0, 1)){
  n     <- length(obs)
  obs01 <- (obs - xlim[1]) / (xlim[2] - xlim[1])

  xx    <- outer(obs01, obs01, FUN = "-")
  
  Kast <- Kast(xx / h)
  a <- 1 / h / n^2 * sum(Kast)
  b <- 2 / n / h * dnorm(0)
  return(a + b) 
} 

# riskKDE2 <- function(obs, xim = c(0, 1), ngrid = 2^15){
#   n <- length(obs)
#   
#   obs01  <- (obs - xlim[1]) / (xlim[2] - xlim[1])
#   
#   fitall <- density(obs01, bw = h, n = ngrid)
# 
#   int2 <- sum(fitall$y^2) * diff(fitall$x)[1]
#   
#   lapply(1:n, function(i) {
#     Xi <- obs01[i]
#     density(obs01[-1])
#     
#   })
#   sumi <- 
#   
#   return(int2 - 2 / n * sumi)
# }

hgrid <- seq(.01, 3, length.out = 30)
J2    <- numeric(length(hgrid))
for (h in seq_along(hgrid)) J2[h] <- riskKDE(Observations, hgrid[h])
layout(1)
plot(hgrid, J2, type = 'l')









# CV_kern <- function(x) {
#   n <- length(x)
#   a <- min(x)
#   b <- max(x)
#   b <- b + (b - a) / n
#   h_CV <- (b - a) / n
# 
#   J0 = 20*dnorm(0)/(n*h_CV)
#   #xx=((1:n)*0+1)%*%t((x))-t(((1:n)*0+1)%*%t((x)))
#   xx = -outer(x, x, FUN = "-")
# 
#   N <- round(ifelse(n < 100, n / 2, n / 4))
# #  if (n<100)
# #	N=round(n/2)
# #else
# #	N=round(n/4)
# #end;
# 
#   J = 1:N;
# 
#   for (m in 1:N) {
#     h = m * (b - a)/N
#     J[m] = 2 * dnorm(0) / (n * h) + (1 / (n^2 * h)) *
#            sum(dnorm(xx / h, 0, sqrt(2)) - 2 * dnorm(xx / h))
#     if (J[m] < J0) {
#       h_CV = h
#       J0 = J[m]
#       }
#   }
# op <- par(mfcol = c(1, 2), pty = "m", omi = c(0, 0, 0, 0))
# 
# plot((1:N)*(b - a) / N, J, type = 'l', lwd = 2, col = 'darkred',
#      main = 'La courbe de la fonction de validation croisée',
#      xlab = 'fenêtre', ylab = 'CV')
# 
# n = length(x)
# a = min(x)
# b = max(x)
# a = a - (b - a)
# b = b + (b - a)
# tt = (a + (b - a)*(1:500)/500) %*% t ((1:n) * 0 + 1)
# xx = ((1:500) * 0 + 1) %*% t(x)
# z = (xx - tt) / h_CV
# hatf = (1 / (n * h_CV)) * as.vector(rowSums(dnorm(z)))
# plot(a + (b - a) * (1:500) / 500, hatf, type = 'l', col = 'darkred',
#      lwd = 2, main = 'Estimateur à noyau avec la fenêtre optimale',
#      xlab = "", ylab = "")
# par(op)
# 
# return(2.2*h_CV)
# }
# 
# 
# CV_kern(Observations)

