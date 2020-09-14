setwd("/Users/august/Desktop/Dissertation")
complete_data <- read.csv("idl_complete_clean.csv")
data <- subset(complete_data, VALIDATION %in% c("Yes", "YES"))
data$DDATE<-as.Date(data$DDATE, "%d/%m/%Y")
data$BDATE<-as.Date(data$BDATE, "%d/%m/%Y")
data <- data[, 2:11]
data$BREGION <- NULL
data$DREGION <- NULL
data$total_age <- data$AGEYEARS + data$DAYSSINCEBD/365
UK <- subset(data, DCOUNTRY == "EW")

bootstrap_est <- function(region, threshold, B){
  
  year_1 <- floor(threshold)
  month_1 <- floor((threshold - year_1) * 12)
  day_1 <- floor(((threshold - year_1) * 12 - month_1) * 30)
  region$achieve_threshold<- region$BDATE %m+% years(year_1)
  region$achieve_threshold <- region$achieve_threshold %m+% months(month_1)
  region$achieve_threshold <- region$achieve_threshold %m+% days(day_1) 
  region$excess_life <- region$total_age - threshold
  region <- subset(region, total_age >= threshold)
  
  n <- length(region$total_age)
  gamma <- rep(0, B)
  sigma_GPD <- rep(0, B)
  sigma_EXP <- rep(0, B)
  p_val <- rep(0, B)
  
  for (i in 1:B){
    index <- sample(1:n, n, replace = TRUE)
    new_region <- region[index, ]
    opt0 <- optim(c(0.5, 0.1), LLGPD_Unbiased, region = new_region)
    opt1 <- optim(0, LLEXP_Unbiased, region = new_region, lower = 0.1, upper = 10, method = "Brent")
    gamma[i] <- opt0$par[2]
    sigma_GPD[i] <- opt0$par[1]
    sigma_EXP[i] <- opt1$par
    LL0 <- -opt0$value
    LL1 <- -opt1$value
    LR <- -2 * (LL1 - LL0)
    p_val[i] <- 1 - pchisq(LR, df = 1)
  }
  
  return (list(gamma = gamma, sigma_GPD = sigma_GPD, sigma_EXP = sigma_EXP, p_val = p_val))

}



p_val_boot <- function(region){
  
  thresholds <- seq(105, 113, 0.1)
  n <- length(thresholds)
  p_val_MLE <- rep(0, n)
  for (i in 1:n){
    print (i)
    threshold <- thresholds[i]
    est <- bootstrap_est(region, threshold, 10)
    p_vals <- est$p_val
    p_val_MLE[i] <- mean(p_vals)
    
    
  }
  
  return (p_val_MLE)
}




LLGPD <- function(paras, region){
  
  sigma <- paras[1]
  gamma <- paras[2]
  x <- region$excess_life
  n <- length(x)
  eps = 1e-9
  
  if (abs(gamma) < eps){
    LL = -n*log(sigma)-sum(x)/sigma
  }
  else{
    temp <- pmax(0, 1+(gamma/sigma)*x)
    LL = -n*log(sigma)-(1/gamma + 1)*sum(log(temp))
  }
  
  return (LL)
}

LLEXP <- function(sigma, region){
  
  x <- region$excess_life
  n <- length(x)
  LL <- -n * log(sigma) - sum(x)/sigma
  
  return (LL)
  
}

LLGPD_Unbiased <- function(paras, region){
  
  sigma <- paras[1]
  gamma <- paras[2]
  
  time <- sort(region$DDATE)
  n <- length(time)
  start_date <- time[1]
  month(start_date) <- 1
  day(start_date) <- 1
  end_date <- time[n]
  year(end_date) <- year(end_date) + 1
  month(end_date) <- 1
  day(end_date) <- 1
  
  group_A <- subset(region, achieve_threshold < start_date)
  group_B <- subset(region, achieve_threshold >= start_date)
  
  r_A <- as.numeric((end_date - group_A$achieve_threshold)/365)
  r_B <- as.numeric((end_date - group_B$achieve_threshold)/365)
  l_A <- as.numeric((start_date - group_A$achieve_threshold)/365)
  n_A <- length(r_A)
  n_B <- length(r_B)
  
  if (n_A == 0 & n_B != 0){
    
    LL <- LLGPD(paras, region) - sum(log(GPD_cdf(paras, r_B)))
    
  }else if(n_A != 0 & n_B == 0){
    
    LL <- LLGPD(paras, region) - sum(log(GPD_cdf(paras, r_A) - GPD_cdf(paras, l_A)))
    
    
  }else{
    
    LL <- LLGPD(paras, region) - sum(log(GPD_cdf(paras, r_B))) - sum(log(GPD_cdf(paras, r_A) - GPD_cdf(paras, l_A)))
    
  }
  
  return (-LL)
}

LLEXP_Unbiased <- function(paras, region){
  
  
  time <- sort(region$DDATE)
  n <- length(time)
  start_date <- time[1]
  month(start_date) <- 1
  day(start_date) <- 1
  end_date <- time[n]
  year(end_date) <- year(end_date) + 1
  month(end_date) <- 1
  day(end_date) <- 1
  
  group_A <- subset(region, achieve_threshold < start_date)
  group_B <- subset(region, achieve_threshold >= start_date)
  
  r_A <- as.numeric((end_date - group_A$achieve_threshold)/365)
  r_B <- as.numeric((end_date - group_B$achieve_threshold)/365)
  l_A <- as.numeric((start_date - group_A$achieve_threshold)/365)
  n_A <- length(r_A)
  n_B <- length(r_B)
  
  
  if (n_A == 0 & n_B != 0){
    
    LL <- LLEXP(paras, region) - sum(log(EXP_cdf(paras, r_B)))
    
  }else if(n_A != 0 & n_B == 0){
    
    LL <- LLEXP(paras, region) - sum(log(EXP_cdf(paras, r_A) - EXP_cdf(paras, l_A)))
    
  }else{
    
    LL <- LLEXP(paras, region) - sum(log(EXP_cdf(paras, r_B))) - sum(log(EXP_cdf(paras, r_A) - EXP_cdf(paras, l_A)))
    
  }
  
  
  return (-LL)
}

GPD_cdf <- function(paras, x){
  
  sigma <- paras[1]
  gamma <- paras[2]
  n <- length(x)
  cdf <- rep(0, n)
  eps = 1e-9
  if (gamma < 0){
    upper <- -sigma/gamma
  }else{
    upper <- Inf
  }
  cdf[x >= upper] <- 1
  x_1 <- x[x < upper]
  
  if (abs(gamma) < eps){
    cdf_1 <- 1 - exp(-x_1/sigma)
  }
  else{
    cdf_1 <- 1 - (1 + gamma * x_1/sigma)^(-1/gamma)
  }
  cdf[x < upper] <- cdf_1
  return (cdf)
}

EXP_cdf <- function(sigma, x){
  
  cdf <- 1 - exp(-x/sigma)
  
  return (cdf)
}

p_value_dis <- function(region, threshold){
  
  year_1 <- floor(threshold)
  month_1 <- floor((threshold - year_1) * 12)
  day_1 <- floor(((threshold - year_1) * 12 - month_1) * 30)
  region$achieve_threshold<- region$BDATE %m+% years(year_1)
  region$achieve_threshold <- region$achieve_threshold %m+% months(month_1)
  region$achieve_threshold <- region$achieve_threshold %m+% days(day_1) 
  region$excess_life <- region$total_age - threshold
  region <- subset(region, total_age >= threshold)
  opt0 <- optim(c(1, 0.1), LLGPD_Unbiased, region = region)
  opt1 <- optim(1, LLEXP_Unbiased, region = region, lower = 0.1, upper = 10, method = "Brent")
  LL0 <- -opt0$value
  LL1 <- -opt1$value
  LR <- -2 * (LL1 - LL0)
  p <- 1 - pchisq(LR, df = 1)
  
  return (p)
  
}
