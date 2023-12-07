rm(list = ls())

#options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("doSNOW",lib = "~/code")

#library()
#.libPaths()
#.libPaths("~/code")

#setwd("~/code")

library("doSNOW")
library("midasr")
library("foreach")
library("doParallel")
library("parallel")

set.seed(100)

lnl_l <- function(paras, y, z, x) {
  the <- paras[1:2]
  wet <- nealmon(p = c(1, the), d = 3)
  xx <- mls(x, 0:2, 3)
  xx0 <- as.matrix(xx) %*% wet
  con <- rep(1, length(y))
  data <- cbind(con, z, xx0)
  beta <- paras[3:5]
  eta <- as.numeric(data %*% beta)
  mll <- sum(y * eta - log(1 + exp(eta)))- sum(beta^2)/8
  return(-mll)
}

lnl <- function(paras, y, z, x, q, ga) {
  the <- paras[1:2]
  neg.part <- ifelse(q <= ga, 1, 0)
  pos.part <- ifelse(q > ga, 1, 0)
  con1 <- neg.part * rep(1, length(y))
  con2 <- pos.part * rep(1, length(y))
  z1 <- neg.part * z
  z2 <- pos.part * z
  xx <- mls(x, 0:2, 3)
  x1 <- neg.part * xx
  x2 <- pos.part * xx
  wet <- nealmon(p = c(1, the), d = 3)
  x10 <- as.matrix(x1) %*% wet
  x20 <- as.matrix(x2) %*% wet
  data <- cbind(con1, z1, x10, con2, z2, x20)
  beta <- paras[3:8]
  eta <- as.numeric(data %*% beta)
  mll <- sum(y * eta - log(1 + exp(eta)))- sum(beta^2)/8
  return(-mll)
}




tmidas_logit <- function(y,x,z,q,n) {
  # linear model estimates
  paras0_l <- c(0,0,0,0,0)
  fit_l = tryCatch(
    {
      nlm(lnl_l, paras0_l, y, z, x)},
    #   # },warning = function(w) {
    #   #   minimum <- Inf
    #   #   estimate <- rep(NA, 8)
    #   #   iterations <- 0
    #   #   return(list(minimum=minimum, estimate=estimate, iterations=iterations))
    #   # },
    error = function(e) {
      minimum <- Inf
      estimate <- rep(NA, 5)
      iterations <- 0
      return(list(minimum=minimum, estimate=estimate, iterations=iterations))
    }
  )
  behat_l = fit_l$estimate[3:5]
  thehat_l = fit_l$estimate[1:2]
  mll_l = -fit_l$minimum
  
  wet <- nealmon(p = c(1, thehat_l), d = 3)
  xx <- mls(x, 0:2, 3)
  xx0 <- as.matrix(xx) %*% wet
  con <- rep(1, length(y))
  data <- cbind(con, z, xx0)
  eta <- as.numeric(data %*% behat_l)
  fitval = 1 / (1 + exp(-eta))
  
  # threshold model estimates
  r1 <- quantile(q, probs = 0.15)
  r2 <- quantile(q, probs = 0.85)
  gas <- seq(r1, r2, by = 0.1)
  lenga = length(gas)
  mlls <- matrix(NA, lenga, 1)
  betas <- matrix(NA, lenga, 6)
  thetas <- matrix(NA, lenga, 2)
  for (j in 1:lenga) {
    #j = 10
    ga <- gas[j]
    paras0 <- c(0,0,0,0,0,0,0,0)
    fit <- tryCatch(
      {
        nlm(lnl, paras0, y, z, x, q, ga)},
      #   # },warning = function(w) {
      #   #   minimum <- Inf
      #   #   estimate <- rep(NA, 8)
      #   #   iterations <- 0
      #   #   return(list(minimum=minimum, estimate=estimate, iterations=iterations))
      #   # },
      error = function(e) {
        minimum <- Inf
        estimate <- rep(NA, 8)
        iterations <- 0
        return(list(minimum=minimum, estimate=estimate, iterations=iterations))
      }
    )
    behat = fit$estimate[3:8]
    thehat = fit$estimate[1:2]
    lnlmin = fit$minimum
    
    mlls[j,] <- lnlmin
    thetas[j, ] <- thehat
    betas[j, ] <- behat
  }
  
  gi <- which.min(mlls)
  gahat <- gas[gi]
  lnlmin <- mlls[gi]
  mllmax = -lnlmin
  behat <- betas[gi, ]
  thehat <- thetas[gi, ]
  
  return(list(behat_l = behat_l, thehat_l = thehat_l, mll_l = mll_l, fitval = fitval, behat = behat, thehat = thehat, gahat = gahat, mllmax = mllmax))
}

tmidas_test1 = function(y,x,z,q,n) {
  tm = tmidas_logit(y,x,z,q,n)
  mll_linear = tm$mll_l
  mll_tmidas = tm$mllmax
  lr1 = 2 * (mll_tmidas - mll_linear)
  
  p = tm$fitval
  
  nrep = 100
  te1 <- matrix(rep(NA, 1 * nrep), ncol = 1)
  
  for (i in 1:nrep) {
    ru = runif(n)
    yb = ifelse(ru < p, 1, 0)
    tmb = tmidas_logit(yb,x,z,q,n)
    mll_linearb = tmb$mll_l
    mll_tmidasb = tmb$mllmax
    lrb = 2 * (mll_tmidasb - mll_linearb)
    te1[i,]=lrb
    
  }
  pvalue = colMeans(te1 >= lr1)
  return(list(lr1=lr1, pvalue=pvalue))
}

# time1 <- Sys.time()
# n <- 1000
# x <- rnorm(3 * n)
# q <- rnorm(n)
# z <- rnorm(n)
# fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
# wx <- 1 + 1 * z + 1 * mls(x, 0:2, 3) %*% fn_x +
#   1 * (q > 0.5) + 1 * (q > 0.5) * z + 1 * (q > 0.5) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
# y <- ifelse(wx > 0, 1, 0)
# tm = tmidas_logit(y,x,z,q,n)
# time2 <- Sys.time()
# print(time2- time1)


srep <- 1000
simp1 <- matrix(rep(NA, 1 * srep), ncol = 1)

cl <- makeSOCKcluster(39)
#clusterEvalQ(cl, .libPaths("~/code"))
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=srep, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# n=100, size
time1 <- Sys.time()
result1 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                   .combine='rbind', .errorhandling = "pass") %dopar% {
                     n <- 100
                     x <- rnorm(3 * n)
                     q <- rnorm(n)
                     z <- rnorm(n)
                     fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                     wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n)
                     y <- ifelse(wx > 0, 1, 0)
                     test1 <- tmidas_test1(y,x,z,q,n)
                     test1$pvalue
                     
                   }
time2 <- Sys.time()
print(time2- time1)
close(pb)

size1 = colMeans(result1 < 0.05)
size1
hist(result1)
hist(result1, breaks = seq(0,1,0.01))


# n=100, power of small threshold effect
time1 <- Sys.time()
results1 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 100
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        0.5 * (q > 0.1) + 0.5 * (q > 0.1) * z + 0.5 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powers1 = colMeans(results1 < 0.05)
powers1
hist(results1, breaks = seq(0,1,0.01))

# n=100, power of medium threshold effect
time1 <- Sys.time()
resultm1 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 100
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        1 * (q > 0.1) + 1 * (q > 0.1) * z + 1 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerm1 = colMeans(resultm1 < 0.05)
powerm1
hist(resultm1, breaks = seq(0,1,0.01))

# n=100, power of large threshold effect
time1 <- Sys.time()
resultl1 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 100
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        2 * (q > 0.1) + 2 * (q > 0.1) * z + 2 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerl1 = colMeans(resultl1 < 0.05)
powerl1
hist(resultl1, breaks = seq(0,1,0.01))

# n=200, size
time1 <- Sys.time()
result2 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                   .combine='rbind', .errorhandling = "pass") %dopar% {
                     n <- 200
                     x <- rnorm(3 * n)
                     q <- rnorm(n)
                     z <- rnorm(n)
                     fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                     wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n)
                     y <- ifelse(wx > 0, 1, 0)
                     test1 <- tmidas_test1(y,x,z,q,n)
                     test1$pvalue
                     
                   }
time2 <- Sys.time()
print(time2- time1)
close(pb)

size2 = colMeans(result2 < 0.05)
size2
hist(result2)
hist(result2, breaks = seq(0,1,0.01))


# n=200, power of small threshold effect
time1 <- Sys.time()
results2 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 200
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        0.5 * (q > 0.1) + 0.5 * (q > 0.1) * z + 0.5 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powers2 = colMeans(results2 < 0.05)
powers2
hist(results2, breaks = seq(0,1,0.01))

# n=200, power of medium threshold effect
time1 <- Sys.time()
resultm2 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 200
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        1 * (q > 0.1) + 1 * (q > 0.1) * z + 1 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerm2 = colMeans(resultm2 < 0.05)
powerm2
hist(resultm2, breaks = seq(0,1,0.01))

# n=200, power of large threshold effect
time1 <- Sys.time()
resultl2 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 200
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        2 * (q > 0.1) + 2 * (q > 0.1) * z + 2 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerl2 = colMeans(resultl2 < 0.05)
powerl2
hist(resultl2, breaks = seq(0,1,0.01))


# n=500, size
time1 <- Sys.time()
result5 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                   .combine='rbind', .errorhandling = "pass") %dopar% {
                     n <- 500
                     x <- rnorm(3 * n)
                     q <- rnorm(n)
                     z <- rnorm(n)
                     fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                     wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n)
                     y <- ifelse(wx > 0, 1, 0)
                     test1 <- tmidas_test1(y,x,z,q,n)
                     test1$pvalue
                     
                   }
time2 <- Sys.time()
print(time2- time1)
close(pb)

size5 = colMeans(result5 < 0.05)
size5
hist(result5)
hist(result5, breaks = seq(0,1,0.01))


# n=500, power of small threshold effect
time1 <- Sys.time()
results5 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 500
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        0.5 * (q > 0.1) + 0.5 * (q > 0.1) * z + 0.5 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powers5 = colMeans(results5 < 0.05)
powers5
hist(results5, breaks = seq(0,1,0.01))

# n=500, power of medium threshold effect
time1 <- Sys.time()
resultm5 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 500
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        1 * (q > 0.1) + 1 * (q > 0.1) * z + 1 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerm5 = colMeans(resultm5 < 0.05)
powerm5
hist(resultm5, breaks = seq(0,1,0.01))

# n=500, power of large threshold effect
time1 <- Sys.time()
resultl5 <- foreach(i=1:srep, .packages="midasr", .options.snow=opts,
                    .combine='rbind', .errorhandling = "pass") %dopar% {
                      n <- 500
                      x <- rnorm(3 * n)
                      q <- rnorm(n)
                      z <- rnorm(n)
                      fn_x <- nealmon(p = c(1, -0.5, 0), d = 3)
                      wx <- -1 - 1 * z - 1 * mls(x, 0:2, 3) %*% fn_x +
                        2 * (q > 0.1) + 2 * (q > 0.1) * z + 2 * (q > 0.1) * mls(x, 0:2, 3) %*% fn_x + 1 * rlogis(n) # 修改
                      y <- ifelse(wx > 0, 1, 0)
                      test1 <- tmidas_test1(y,x,z,q,n)
                      test1$pvalue
                      
                    }
time2 <- Sys.time()
print(time2- time1)
close(pb)

powerl5 = colMeans(resultl5 < 0.05)
powerl5
hist(resultl5, breaks = seq(0,1,0.01))





stopCluster(cl)
