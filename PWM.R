PWM_solver <- function(region, k){
  
  x <- region$excess_life
  u_1 <- 1
  v_1 <- 0
  u_2 <- 0
  v_2 <- 1
  m1 <- pwMoment(x, j = u_1, k = v_1)
  m2 <- pwMoment(x, j = u_2, k = v_2)
  
  
  a <- beta(u_1 + 1, v_1 + 1) - beta(u_1 + 1, v_1 + k + 1)
  b <- beta(u_2 + 1, v_2 + 1) - beta(u_2 + 1, v_2 + k + 1)
  
  y <- a/b - m1/m2
  
  return (y)
  
  
}

PWM_est <- function(region){
  
  u_1 <- 1
  v_1 <- 0
  u_2 <- 0
  v_2 <- 1
  x <- region$excess_life
  m1 <- pwMoment(x, j = u_1, k = v_1)
  k <- nleqslv(0.4, PWM_solver, region = region)[[1]]
  sigma <- k * m1/(beta(u_1 + 1, v_1 + 1) - beta(u_1 + 1, v_1 + 1 + k))
  
  return (c(sigma, -k))
  
  
}



plot_hist <-function(){
  
  par(mfrow=c(1,2))
  hist(UK$excess_life, xlab = "Excess Age (Years)", ylab = "Count")
  hist(FRA$excess_life, xlab = "Excess Age (Years)", ylab = "Count")
    
  

}

