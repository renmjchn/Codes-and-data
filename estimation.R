
rm(list = ls())

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# install.packages("mcmc",lib = "~/code")

# library()
# .libPaths()
#.libPaths("~/code")

# setwd("~/code")

library("doSNOW")
library("midasr")
library("foreach")
library("doParallel")
library("parallel")

set.seed(100)

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

logit <- function() {
  ########################################################
  n <- 500
  x <- rnorm(3 * n)
  q <- rnorm(n)
  z <- rnorm(n)
  ##############################################################
  fn_x <- nealmon(p = c(1, -5, -5), d = 3)
  wx <- -1 - 1*z - 1*mls(x, 0:2, 3)%*%fn_x +
    2*(q > 0.1) + 2*(q > 0.1)*z + 2*(q > 0.1)*mls(x, 0:2, 3)%*%fn_x + rlogis(n) # 修改
  y <- ifelse(wx > 0, 1, 0)
  # 估计程序
  r1 <- quantile(q, probs = 0.15)
  r2 <- quantile(q, probs = 0.85)
  gas <- seq(r1, r2, by = 0.1)
  lenga = length(gas)
  mlls <- matrix(NA, lenga, 1)
  betas <- matrix(NA, lenga, 6)
  thetas <- matrix(NA, lenga, 2)
  wets = matrix(NA, lenga, 3)
  for (j in 1:lenga) {
    #j = 10
    ga <- gas[j]
    paras0 <- c(0,0,0,0,0,0,0,0)
    #fit = nlm(lnl, paras0, y, z, x, q, ga)
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
    lnlmin <- fit$minimum
    mlls[j,] <- lnlmin
    thetas[j, ] <- thehat
    betas[j, ] <- behat
    wets[j, ] <- nealmon(p = c(1, thehat[1], thehat[2]), d = 3)
  }
  
  gi <- which.min(mlls)
  gi
  gahat <- gas[gi]
  gahat
  lnlmin <- mlls[gi]
  behat <- betas[gi, ]
  behat
  thehat <- thetas[gi, ]
  thehat
  wethat = wets[gi, ]
  wethat
  plot(mlls)
  return(c(gahat, wethat, thehat, behat))
}



srep <- 1000

cl <- makeSOCKcluster(39)
#clusterEvalQ(cl, .libPaths("~/code"))
registerDoSNOW(cl)

pb <- txtProgressBar(min = 1, max = srep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time1 <- Sys.time()
result <- foreach(
  i = 1:srep, .packages = c("midasr"), .options.snow = opts,
  .combine = "rbind", .errorhandling = "pass"
) %dopar% {
  logit()
}
time2 <- Sys.time()
print(time2 - time1)

close(pb)
stopCluster(cl)

# result = result[-355,]
# result = as.matrix(result)
# result = result[result[,1]<1000,]
meanre <- colMeans(result)
round(meanre, 3)
sdre <- cbind(sd(result[, 1]), sd(result[, 2]), sd(result[, 3]), sd(result[, 4]), sd(result[, 5]), sd(result[, 6]), sd(result[, 7]), sd(result[, 8]), sd(result[, 9]), sd(result[, 10]), sd(result[, 11]), sd(result[, 12]))
round(sdre, 3)
hist(result[,1],main = "gamma")
hist(result[,2],main = "wet1")
hist(result[,3],main = "wet2")
hist(result[,4],main = "wet3")
hist(result[,5],main = "the1")
hist(result[,6],main = "the2")
hist(result[,7],main = "b10")
hist(result[,8],main = "b11")
hist(result[,9],main = "b12")
hist(result[,10],main = "b20")
hist(result[,11],main = "b21")
hist(result[,12],main = "b22")


