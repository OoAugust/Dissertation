MH <- function(region, burn_in, iter){
  
  N <- burn_in + iter
  x <- region$excess_life
  paras <- c(0, 0.1)
  output <- list(1, 0.1)
  
  for (i in 1:N){
    
    
    eta <- paras[1]
    gamma <- paras[2]
    eps_gamma <- rnorm(1, mean= 0, sd = 0.1)
    eps_eta <- rnorm(1, mean = 0, sd = 1)
    new_gamma <- gamma + eps_gamma
    new_eta <- eta + eps_eta
    new_para <- c(new_eta, new_gamma)
    
    log_acc_prob <- min(0, log(prior(new_para)) + LL(new_para, x) - log((prior(paras))) - LL(paras, x))
    

    u <- runif(1, min = 0, max = 1)
    
    if (log(u) <= log_acc_prob & i > burn_in & !is.na(log_acc_prob)){
      output[[1]] <- c(output[[1]], exp(new_eta))
      output[[2]] <- c(output[[2]], new_gamma)
      paras <- new_para
    }
    
  }
  
  return (output)
}

prior <- function(paras){
  
  p <- dnorm(paras[1], 0, 10) * dnorm(paras[2], 0, 1)
  
  return (p)
  
}

LL <- function(paras, x){
  
  eta <- paras[1]
  gamma <- paras[2]
  sigma <- exp(eta)
  eps = 1e-9
  n <- length(x)
  
  if (abs(gamma) < eps){
    LL = -n*log(sigma)-sum(x)/sigma
  }
  else{
    LL = -n*log(sigma)-(1/gamma + 1)*sum(log(pmax(0, 1+(gamma/sigma)*x)))
  }
  
  return (LL)
  
  
}

plot_bay <- function(region){
  
  bayesian<- MH(region, 10000, 100000)
  print (c(mean(bayesian[[1]]), var(bayesian[[1]])))
  print (c(mean(bayesian[[2]]), var(bayesian[[2]])))
  sigma <- bayesian[[1]]
  xi <- bayesian[[2]]
  n <- length(sigma)
  par(mfrow=c(1,2))
  plot(c(1:n), sigma, xlim = c(0, 400), "l", ylab = expression(sigma), xlab = "Count", col = "skyblue", lty = 1)
  abline(h = mean(sigma), col = "red", lwd=1, lty=1)
  legend("bottomright", col = c("skyblue", "red"), legend = c("Posterior Samples", "Posterior Mean"), lty = c(1, 1))
  plot(c(1:n), xi, xlim = c(0, 400), "l", ylab = expression(xi), xlab = "Count", col = "skyblue", lty = 1)
  abline(h = mean(xi), col = "red", lwd=1, lty=1)
  legend("topright", col = c("skyblue", "red"), legend = c("Posterior Samples", "Posterior Mean"), lty = c(1, 1))
}

new_est <- function(){
  
  thresholds <- seq(105, 112, 0.1)
  n <- length(thresholds)
  sigma_bay <- rep(0, n)
  xi_bay <- rep(0, n)
  sigma_MLE <- rep(0, n)
  xi_MLE <- rep(0, n)
  for (i in 1:n){
    
    print (i)
    threshold <- thresholds[i]
    world <- set_threshold(threshold)
    NN <- subset(world, DCOUNTRY %in% c("DNK", "DEU", "BEL", "FIN", "NOR", "SWE", "CHE"))
    bayesian <- MH(NN, 10000, 20000)
    sigma_bay[i] <- mean(bayesian[[1]])
    xi_bay[i] <- mean(bayesian[[2]])
    
    est <- optim(c(1, 0.1), LLGPD, region = NN)
    sigma_MLE[i] <- est$par[[1]]
    xi_MLE[i] <- est$par[[2]]
    
  }
  
  return (list(sigma_bay = sigma_bay, xi_bay = xi_bay, sigma_MLE = sigma_MLE, xi_MLE = xi_MLE))
  
  
}

plot_new <- function(){
  
  ets <- new_est()
  thresholds <- seq(105, 112, 0.1)
  
  par(mfrow=c(1,2))
  
  plot(thresholds, a[[1]], xlim = c(105, 112), "l", ylab = expression(sigma), xlab = "Age", col = "skyblue", lty = 1, ylim = c(0.5, 2))
  par(new = TRUE)
  plot(thresholds, a[[3]], xlim = c(105, 112), "l", ylab = expression(sigma), xlab = "Age", col = "red", lty = 1, ylim = c(0.5, 2))
  legend("topleft", col = c("skyblue", "red"), legend = c("Bayesian Estimates", "MLE Estimates"), lty = c(1, 1), cex =  1.3)
  
  plot(thresholds, a[[2]], xlim = c(105, 112), "l", ylab = expression(xi), xlab = "Age", col = "skyblue", lty = 1, ylim = c(-1.5, 0.5))
  par(new = TRUE)
  plot(thresholds, a[[4]], xlim = c(105, 112), "l", ylab = expression(xi), xlab = "Age", col = "red", lty = 1, ylim = c(-1.5, 0.5))
  legend("topleft", col = c("skyblue", "red"), legend = c("Bayesian Estimates", "MLE Estimates"), lty = c(1, 1), cex = 1.3)
  
  
  
}