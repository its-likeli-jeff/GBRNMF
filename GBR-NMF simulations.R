#regular nmf simulation
set.seed(62341)
simh <- matrix(runif(800,3,8), 2, 400)
simw <- matrix(runif(2000,3,8), 1000, 2)
simx <- simw %*% simh + matrix(rnorm(400*1000),1000,400)
set.seed(5108651)
simrun50k <- nmf(simx, 2, maxit=500000, eps=0)


head(simw)
head(simrun50k$w)

head(t(simh))
head(t(simrun50k$h))

#gbr_nmf simulation
set.seed(576)

total_w_cor <- matrix(0, nrow = 10, ncol = 7)
total_a_cor <- vector("numeric", length= 7)
total_s_cor <- matrix(0, nrow = 10, ncol = 7)

for(k in 1:10){
  sims <- matrix(runif(14000, 0, 1), 7, 2000)

  simw_1 <- cbind(rnorm(100,mean = 1, sd=0.1),matrix(0,nrow = 100,ncol = 3), rnorm(100,mean = 1, sd=0.1), rnorm(100,mean = 1, sd=0.1), rnorm(100,mean = 1, sd=0.1))
  simw_2 <- cbind(matrix(0, nrow = 100, ncol = 1), rnorm(100, mean = 1, sd = 0.25), matrix(0, nrow = 100, ncol = 2), rnorm(100,mean = 1, sd=0.05), rnorm(100,mean = 1, sd=0.3), rnorm(100, mean = 1, sd = 0.05))
  simw_3 <- cbind(matrix(0, nrow = 100, ncol = 2), rnorm(100, mean = 1, sd = 0.15), matrix(0, nrow = 100, ncol = 1), rnorm(100,mean = 1, sd=0.05), rnorm(100,mean = 1, sd=0.2), rnorm(100, mean = 1, sd = 0.25))
  simw_4 <- cbind(matrix(0, nrow = 100, ncol = 3), rnorm(100, mean = 1, sd = 0.05), rnorm(100,mean = 1, sd=0.2), rnorm(100,mean = 1, sd=0.15), rnorm(100, mean = 1, sd = 0.3))
  simw <- rbind(simw_1,simw_2,simw_3,simw_4)
  
  #make a shift for W. Explain in paper. Include a diagram of some sorts
  #generate columns of W binary as well possibly
  #can do both
  #add a diagram

  sima <- diag(rnorm(7,mean = 2, sd = 1.5))

  simx <- simw %*% sima %*% sims
  
  groups <- rep(1:4,each=100)
  gmat <- matrix(0, nrow = 400, ncol = 4)
  for(i in 1:400){
    gmat[i,groups[i]] <- 1
  }

  constrainedrun <- gbr_nmf(x = simx, q = 7, maxit = 50000, w = gmat, s = t(sims[5,]))

  #Correlation for W
  for(i in 1:ncol(simw)){
    print(cor(simw[,i], constrainedrun$w[,i]))
    total_w_cor[k,i] <- cor(simw[,i], constrainedrun$w[,i])
  }

  #Correlation for A
  cor(as.vector(constrainedrun$a), as.vector(sima))
  total_a_cor[k] <- cor(as.vector(constrainedrun$a), as.vector(sima))
  
  #Correlations for S
  for(i in 1:(nrow(sims))){
    print(cor(sims[i,],constrainedrun$s[i,]))
    total_s_cor[k,i] <- cor(sims[i,],constrainedrun$s[i,])
  }
}

min(simw)
max(simw)
