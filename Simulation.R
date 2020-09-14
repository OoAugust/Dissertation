diff_xi <- function(){
  
  size <- 500
  xi <- c(-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1)
  n <- length(xi)
  col_names <- as.character(xi)
  n_sim <- 1000 
  
  
  MLE_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_xi) <- col_names
  Bay_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_xi) <- col_names
  PWM_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_xi) <- col_names
  MLE_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_sigma) <- col_names
  Bay_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_sigma) <- col_names
  PWM_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_sigma) <- col_names
  acc_prob <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(acc_prob) <- col_names
  
  sigma <- 1.3
  
  for (i in 1:n_sim){
    print (i)
    for (j in 1:n){
      
      xi_j <- xi[j]
      x <- rgpd(size, scale = sigma, shape = xi_j)
      est_pwm <- PWM_est(x)
      PWM_sigma[i, j] <- est_pwm[1]
      PWM_xi[i, j] <- est_pwm[2]
      
      est_bay <- MH(x, 1000, 5000)
      sigma_bay <- est_bay[[1]]
      xi_bay <- est_bay[[2]]
      N <- length(sigma_bay)
      Bay_sigma[i, j] <- mean(sigma_bay[-1])
      Bay_xi[i, j] <- mean(xi_bay[-1])
      acc_prob[i, j] <- (N - 1)/5000
      
      opt <- optim(c(1, 0.1), LLGPD, x = x)
      MLE_sigma[i, j] <- opt$par[1]
      MLE_xi[i, j] <- opt$par[2]
    }
    
  }
  
  return (list(MLE_xi = MLE_xi, PWM_xi = PWM_xi, Bayesian_xi = Bay_xi, MLE_sigma = MLE_sigma, PWM_sigma = PWM_sigma, Bayesian_sigma = Bay_sigma, acc_prob = acc_prob))
  
  
}


plot_xi <- function(){
  
  sigma <- 1.3
  xi <- c(-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1)
  result2 <- diff_size()
  MLE_xi <- result2[[1]]
  PWM_xi <- result2[[2]]
  Bay_xi <- result2[[3]]
  MLE_sigma <- result2[[4]]
  PWM_sigma <- result2[[5]]
  Bay_sigma <- result2[[6]]
  sizes <- c(10, 25, 50, 70, 100, 200, 300, 500)
  
  par(mfrow=c(1,2))
  plot(xi, colMeans(MLE_xi) - xi, col = "skyblue", lty = 1, xlim = c(-1, 1), ylim = c(-0.2, 0.6), ylab = "Bias", xlab = expression(xi), type = "b", pch = 19)
  lines(xi, colMeans(PWM_xi) - xi, col = "hotpink", lty = 1, xlim =c(-1, 1), ylim = c(-0.2, 0.6), xlab = expression(xi), ylab = "Bias", type = "b", pch = 19)
  lines(xi, colMeans(Bay_xi) - xi, col = "orange", lty = 1, xlim = c(-1, 1), ylim = c(-0.2, 0.6), xlab = expression(xi), ylab = "Bias", type = "b", pch = 19)
  legend("topright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.3)
  
  plot(xi, sqrt(colVars(as.matrix(MLE_xi))), col = "skyblue", lty = 1, xlim = c(-1, 1), ylim = c(0, 0.15), ylab = "Standard Error", xlab = expression(xi), type = "b", pch = 19)
  lines(xi, sqrt(colVars(as.matrix(PWM_xi))), col = "hotpink", lty = 1, xlim = c(-1, 1), ylim = c(0, 0.15), xlab = "Sample Size", ylab = expression(xi), type = "b", pch = 19)
  lines(xi, sqrt(colVars(as.matrix(Bay_xi))), col = "orange", lty = 1, xlim = c(-1, 1), ylim = c(0, 0.15), xlab = "Sample Size", ylab = expression(xi), type = "b", pch = 19)
  
  legend("bottomright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.3)
  
}
  



diff_size <- function(){
  
  sizes <- c(10, 25, 50, 70, 100, 200, 300, 500)
  n <- length(sizes)
  col_names <- as.character(sizes)
  n_sim <- 1000
  
  MLE_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_xi) <- col_names
  Bay_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_xi) <- col_names
  PWM_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_xi) <- col_names
  MLE_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_sigma) <- col_names
  Bay_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_sigma) <- col_names
  PWM_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_sigma) <- col_names
  acc_prob <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(acc_prob) <- col_names
  
  sigma <- 1.3
  xi <- -0.1
  
  for (i in 1:n_sim){
    print (i)
    for (j in 1:n){
      
      size <- sizes[j]
      x <- rgpd(size, scale = sigma, shape = xi)
      est_pwm <- PWM_est(x)
      PWM_sigma[i, j] <- est_pwm[1]
      PWM_xi[i, j] <- est_pwm[2]
      
      est_bay <- MH(x, 1000, 5000)
      sigma_bay <- est_bay[[1]]
      xi_bay <- est_bay[[2]]
      N <- length(sigma_bay)
      Bay_sigma[i, j] <- mean(sigma_bay[-1])
      Bay_xi[i, j] <- mean(xi_bay[-1])
      acc_prob[i, j] <- (N - 1)/5000
      
      opt <- optim(c(1, 0.1), LLGPD, x = x)
      MLE_sigma[i, j] <- opt$par[1]
      MLE_xi[i, j] <- opt$par[2]
    }
    
  }
  
  return (list(MLE_xi = MLE_xi, PWM_xi = PWM_xi, Bayesian_xi = Bay_xi, MLE_sigma = MLE_sigma, PWM_sigma = PWM_sigma, Bayesian_sigma = Bay_sigma, acc_prob = acc_prob))
  
  
}


diff_threshold <- function(){
  
  
  size <- 10000
  thresholds <- c(102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
  col_names <- as.character(thresholds)
  n <- length(thresholds)
  n_sim <- 500
  
  MLE_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_xi) <- col_names
  Bay_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_xi) <- col_names
  PWM_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_xi) <- col_names
  MLE_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_sigma) <- col_names
  Bay_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_sigma) <- col_names
  PWM_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_sigma) <- col_names
  acc_prob <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(acc_prob) <- col_names
  
  sigma <- 1.3
  xi <- -0.1
  loc <- 105
  
  for (i in 1:n_sim){
    print (i)

    for (j in 1:n){
      
      x <- rgev(size, loc = loc, scale = sigma, shape = xi)
      threshold <- thresholds[j]
      x <- x - threshold
      x <- x[x >= 0]
      if (length(x) <= 1){
        break
      }
      est_pwm <- PWM_est(x)
      PWM_sigma[i, j] <- est_pwm[1]
      PWM_xi[i, j] <- est_pwm[2]
      
      est_bay <- MH(x, 1000, 5000)
      sigma_bay <- est_bay[[1]]
      xi_bay <- est_bay[[2]]
      N <- length(sigma_bay)
      Bay_sigma[i, j] <- mean(sigma_bay[-1])
      Bay_xi[i, j] <- mean(xi_bay[-1])
      acc_prob[i, j] <- (N - 1)/5000
      
      opt <- optim(c(1, 0.1), LLGPD, x = x)
      MLE_sigma[i, j] <- opt$par[1]
      MLE_xi[i, j] <- opt$par[2]
    }
    
  }
  
  return (list(MLE_xi = MLE_xi, PWM_xi = PWM_xi, Bayesian_xi = Bay_xi, MLE_sigma = MLE_sigma, PWM_sigma = PWM_sigma, Bayesian_sigma = Bay_sigma, acc_prob = acc_prob))
  
  
}




plot_threshold <- function(){
  
  sigma <- 1.3
  xi <- -0.1
  result <- diff_size()
  MLE_xi <- result[[1]]
  PWM_xi <- result[[2]]
  Bay_xi <- result[[3]]
  MLE_sigma <- result[[4]]
  PWM_sigma <- result[[5]]
  Bay_sigma <- result[[6]]
  thresholds <- c(103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
  
  sigma_1 <- sigma + xi * (thresholds - 105)
  par(mfrow=c(1,2))
  plot(thresholds, colMeans(MLE_xi, na.rm = TRUE)[-1] - xi, col = "skyblue", lty = 1, xlim = c(103, 112), ylim = c(-1.5, 0.5), ylab = "Bias", xlab = "Threshold", type = "b", pch = 19)
  lines(thresholds, colMeans(PWM_xi, na.rm = TRUE)[-1] - xi, col = "hotpink", lty = 1, xlim =c(103, 112), ylim = c(-1.5, 0.5), xlab = "Threshold", ylab = "Bias", type = "b", pch = 19)
  lines(thresholds, colMeans(Bay_xi, na.rm = TRUE)[-1] - xi, col = "orange", lty = 1, xlim = c(103, 112), ylim = c(-1.5, 0.5), xlab = "Threshold", ylab = "Bias", type = "b", pch = 19)
  legend("topleft", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.2)
  
  plot(thresholds, colMeans(MLE_sigma, na.rm = TRUE)[-1] - sigma_1, col = "skyblue", lty = 1, xlim = c(103, 112), ylim = c(0, 5), ylab = "Bias", xlab = "Threshold", type = "b", pch = 19)
  lines(thresholds, colMeans(PWM_sigma, na.rm = TRUE)[-1] - sigma_1, col = "hotpink", lty = 1, xlim =c(103, 112), ylim = c(0, 5), xlab = "Threshold", ylab = "Bias", type = "b", pch = 19)
  lines(thresholds, colMeans(Bay_sigma, na.rm = TRUE)[-1] - sigma_1, col = "orange", lty = 1, xlim = c(103, 112), ylim = c(0, 5), xlab = "Threshold", ylab = "Bias", type = "b", pch = 19)
  legend("topright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.2)

  
}










plot_size <- function(){
  
  sigma <- 1.3
  xi <- -0.1
  result <- diff_size()
  MLE_xi <- result[[1]]
  PWM_xi <- result[[2]]
  Bay_xi <- result[[3]]
  MLE_sigma <- result[[4]]
  PWM_sigma <- result[[5]]
  Bay_sigma <- result[[6]]
  sizes <- c(10, 25, 50, 70, 100, 200, 300, 500)
  
  par(mfrow=c(1,2))
  plot(sizes, colMeans(MLE_sigma) - sigma, col = "skyblue", lty = 1, xlim = c(0, 500), ylim = c(0, 0.9), ylab = "Bias", xlab = "Sample Size", type = "b", pch = 19)
  lines(sizes, colMeans(PWM_sigma) - sigma, col = "hotpink", lty = 1, xlim =c(0, 500), ylim = c(0, 0.9), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  lines(sizes, colMeans(Bay_sigma) - sigma, col = "orange", lty = 1, xlim = c(0, 500), ylim = c(0, 0.9), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  legend("topright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.3)
  
  plot(sizes, sqrt(colVars(as.matrix(MLE_sigma))), col = "skyblue", lty = 1, xlim = c(0, 500), ylim = c(0, 1.5), ylab = "Standard Error", xlab = "Sample Size", type = "b", pch = 19)
  lines(sizes, sqrt(colVars(as.matrix(PWM_sigma))), col = "hotpink", lty = 1, xlim = c(0, 500), ylim = c(0, 1.5), xlab = "Sample Size", ylab = "Standard Error", type = "b", pch = 19)
  lines(sizes, sqrt(colVars(as.matrix(Bay_sigma))), col = "orange", lty = 1, xlim = c(0, 500), ylim = c(0, 1.5), xlab = "Sample Size", ylab = "Standard Error", type = "b", pch = 19)
  
  legend("topright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.3)
  
}




LLGPD <- function(paras, x){
  
  sigma <- paras[1]
  gamma <- paras[2]
  n <- length(x)
  eps = 1e-9
  
  if (abs(gamma) < eps){
    LL = -n*log(sigma)-sum(x)/sigma
  }
  else{
    temp <- pmax(0, 1+(gamma/sigma)*x)
    LL = -n*log(sigma)-(1/gamma + 1)*sum(log(temp))
  }
  
  return (-LL)
}


PWM_solver <- function(y, k){

  u_1 <- 1
  v_1 <- 0
  u_2 <- 0
  v_2 <- 1
  m1 <- pwMoment(y, j = u_1, k = v_1)
  m2 <- pwMoment(y, j = u_2, k = v_2)
  
  
  a <- beta(u_1 + 1, v_1 + 1) - beta(u_1 + 1, v_1 + k + 1)
  b <- beta(u_2 + 1, v_2 + 1) - beta(u_2 + 1, v_2 + k + 1)
  
  y <- a/b - m1/m2
  
  return (y)
  
  
}

PWM_est <- function(x){
  
  u_1 <- 1
  v_1 <- 0
  u_2 <- 0
  v_2 <- 1
  m1 <- pwMoment(x, j = u_1, k = v_1)
  k <- nleqslv(0.4, PWM_solver, y = x)[[1]]
  sigma <- k * m1/(beta(u_1 + 1, v_1 + 1) - beta(u_1 + 1, v_1 + 1 + k))
  
  return (c(sigma, -k))
  
  
}

MH <- function(x, burn_in, iter){
  
  N <- burn_in + iter
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



diff_threshold_size <- function(){
  
  
  sizes <- c(30, 50, 100, 300, 500, 1000, 3000, 5000)
  threshold <- 108
  col_names <- as.character(sizes)
  n <- length(sizes)
  n_sim <- 1000
  
  MLE_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_xi) <- col_names
  Bay_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_xi) <- col_names
  PWM_xi <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_xi) <- col_names
  MLE_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(MLE_sigma) <- col_names
  Bay_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(Bay_sigma) <- col_names
  PWM_sigma <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(PWM_sigma) <- col_names
  acc_prob <- data.frame(matrix(nrow = n_sim, ncol = n))
  colnames(acc_prob) <- col_names
  
  sigma <- 1.3
  xi <- -0.1
  loc <- 105
  
  for (i in 1:n_sim){
    print (i)
    
    for (j in 1:n){
      
      size<- sizes[j]
      x <- rgev(size, loc = loc, scale = sigma, shape = xi)
      x <- x - threshold
      x <- x[x >= 0]
      if (length(x) <= 1){
        break
      }
      est_pwm <- PWM_est(x)
      PWM_sigma[i, j] <- est_pwm[1]
      PWM_xi[i, j] <- est_pwm[2]
      
      est_bay <- MH(x, 1000, 5000)
      sigma_bay <- est_bay[[1]]
      xi_bay <- est_bay[[2]]
      N <- length(sigma_bay)
      Bay_sigma[i, j] <- mean(sigma_bay[-1])
      Bay_xi[i, j] <- mean(xi_bay[-1])
      acc_prob[i, j] <- (N - 1)/5000
      
      opt <- optim(c(1, 0.1), LLGPD, x = x)
      MLE_sigma[i, j] <- opt$par[1]
      MLE_xi[i, j] <- opt$par[2]
    }
    
  }
  
  return (list(MLE_xi = MLE_xi, PWM_xi = PWM_xi, Bayesian_xi = Bay_xi, MLE_sigma = MLE_sigma, PWM_sigma = PWM_sigma, Bayesian_sigma = Bay_sigma, acc_prob = acc_prob))
  
  
}

plot_threshold_size <- function(){
  
  sigma <- 1.3
  xi <- -0.1
  result <- diff_size()
  MLE_xi <- result[[1]]
  PWM_xi <- result[[2]]
  Bay_xi <- result[[3]]
  MLE_sigma <- result[[4]]
  PWM_sigma <- result[[5]]
  Bay_sigma <- result[[6]]
  threshold <- 108
  sizes <- c(30, 50, 100, 300, 500, 1000, 3000, 5000)
  
  sigma_1 <- sigma + xi * (threshold - 105)
  par(mfrow=c(1,2))
  plot(sizes, colMeans(MLE_xi, na.rm = TRUE) - xi, col = "skyblue", lty = 1, xlim = c(0, 5000), ylim = c(-2, 0.5), ylab = "Bias", xlab = "Sample Size", type = "b", pch = 19)
  lines(sizes, colMeans(PWM_xi, na.rm = TRUE) - xi, col = "hotpink", lty = 1, xlim =c(0, 5000), ylim = c(-2, 0.5), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  lines(sizes, colMeans(Bay_xi, na.rm = TRUE) - xi, col = "orange", lty = 1, xlim = c(0, 5000), ylim = c(-2, 0.5), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  legend("bottomright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.2)
  
  plot(sizes, colMeans(MLE_sigma, na.rm = TRUE) - sigma_1, col = "skyblue", lty = 1, xlim = c(0, 5000), ylim = c(0, 2), ylab = "Bias", xlab = "Sample Size", type = "b", pch = 19)
  lines(sizes, colMeans(PWM_sigma, na.rm = TRUE) - sigma_1, col = "hotpink", lty = 1, xlim = c(0, 5000), ylim = c(0, 2), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  lines(sizes, colMeans(Bay_sigma, na.rm = TRUE) - sigma_1, col = "orange", lty = 1, xlim = c(0, 5000), ylim = c(0, 2), xlab = "Sample Size", ylab = "Bias", type = "b", pch = 19)
  legend("topright", col = c("skyblue", "hotpink", "orange"), legend = c("MLE", "PWM", "Bayesian"), lty = c(1, 1, 1), pch = c(19, 19, 19), cex = 1.2)
  
  
}