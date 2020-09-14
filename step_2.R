uk <- read.csv("uk.csv")
uk$Age <- as.character(uk$Age)
uk$Age[uk$Age == "110+"] <- "111"
uk$Age <- gsub("\\+", "", uk$Age)
uk$Age <- as.numeric(uk$Age)
uk$Year <- paste(uk$Year, "01", sep = "/")
uk$Year <- paste(uk$Year, "01", sep = "/")
uk$Year <- as.Date(uk$Year, "%Y/%m/%d")
uk$birth_year <- uk$Year
year(uk$birth_year) <- year(uk$Year) - uk$Age

#start_date <- "1984-01-01"
#end_date <- "2010-01-01"
#u <- 105  #threshold 
#grouped_data <- process_data(start_date, end_date, u)
#opt_total_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, option = "total", hessian = TRUE)
#opt_male_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, option = "male", hessian = TRUE)
#opt_female_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, option = "female", hessian = TRUE)
#opt_total_EXP <- optim(1, LL_EXP, data = grouped_data, option = "total", hessian = TRUE, method = "Brent", lower = 0, upper = 10)


plot_1 <- function(){
  
  lower <- 90
  upper <- 109#maximum 109, cannot be 109
  ages <- seq(lower, upper, 1)
  n <- length(ages)
  p_values_dis <- rep(0, upper - lower + 1)
  p_values_gender_GPD <- rep(0, upper - lower + 1)
  p_values_gender_EXP <- rep(0, upper - lower + 1)
  start_date <- "1900-01-01"
  end_date <- "2020-01-01"
  
  for (u in lower:upper){
    
    grouped_data <- process_data(start_date, end_date, u)
    
    opt_total_EXP <- optim(1, LL_EXP, data = grouped_data, option = "total", u = u, hessian = TRUE, method = "Brent", lower = 0, upper = 10)
    opt_total_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, u = u, option = "total", hessian = TRUE)
    opt_male_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, u = u, option = "male", hessian = TRUE)
    opt_female_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, u = u, option = "female", hessian = TRUE)
    opt_male_EXP <- optim(1, LL_EXP, data = grouped_data, u = u, option = "male", hessian = TRUE, method = "Brent", lower = 0, upper = 10)
    opt_female_EXP <- optim(1, LL_EXP, data = grouped_data, u = u, option = "female", hessian = TRUE, method = "Brent", lower = 0, upper = 10)
    
    total_LL_GPD <- -opt_total_GPD$value
    total_LL_EXP <- -opt_total_EXP$value
    LL_GPD_male <- -opt_male_GPD$value
    LL_GPD_female <- -opt_female_GPD$value
    LL_EXP_male <- -opt_male_EXP$value
    LL_EXP_female <- -opt_female_EXP$value
    
    LR_dis <- -2 * (total_LL_EXP - total_LL_GPD)
    p_values_dis[u - lower + 1] <- 1 - pchisq(LR_dis, df = 1)
    LR_gender_GPD <- -2 * (total_LL_GPD - LL_GPD_female - LL_GPD_male)
    p_values_gender_GPD[u - lower + 1] <- 1 - pchisq(LR_gender_GPD, df = 2)
    LR_gender_EXP <- -2 * (total_LL_EXP - LL_EXP_female - LL_EXP_male)
    p_values_gender_EXP[u - lower + 1] <- 1 - pchisq(LR_gender_EXP, df = 1)
    
  }
  
  index <- Position(function(x) x > 0.05, p_values_dis)
  p_values_gender <- c(p_values_gender_GPD[1:index - 1], p_values_gender_EXP[index:n])
  
  plot(ages, p_values_dis, col = "orange", type = "b", pch = 19, ylim = c(0, 1), xlim = c(90, 110), ylab = "p values", xlab = "threhold")
  lines(ages, p_values_gender, col = "deepskyblue", type = "b", pch = 19, ylim = c(0, 1), lty = 1, xlim = c(90, 110), ylab = "p values", xlab = "threhold")
  abline(h = 0.05, col = "red", lwd=1, lty=2)
  legend("topleft", col = c("orange", "deepskyblue", "red"), legend = c("Gender", "Distribution", "5% Significance"), lty = c(1, 1, 2))
  
}

plot_2 <- function(u){
  
  start_dates <- seq(as.Date("1811-01-01"), as.Date("1891-01-01"), "20 year")
  end_dates <- seq(as.Date("1831-01-01"), as.Date("1911-01-01"), "20 year")
  n <- length(start_dates)
  sigma_male <- rep(0, n)
  sigma_female <- rep(0, n)
  gamma_male <- rep(0, n)
  gamma_female <- rep(0, n)
  
  for (i in 1:n){
    
    
    start <- start_dates[i]
    end <- end_dates[i]
    grouped_data <- process_data(start, end, u)
    
    opt_male_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, u = u, option = "male")
    opt_female_GPD <- optim(c(1, 0.1), LL_GPD, data = grouped_data, u = u, option = "female")
    est_male <- opt_male_GPD$par
    est_female <- opt_female_GPD$par
    
    sigma_male[i] <- est_male[1]
    gamma_male[i] <- est_male[2]
    sigma_female[i] <- est_female[1]
    gamma_female[i] <- est_female[2]
  
  }
  
  
  par(mfrow=c(1,2))
  
  plot(1:n, smooth(sigma_male), type = "l", pch = 19, ylim = c(1, 3), xlim = c(1, n), col = "blue", ylab = expression(sigma), xlab = "Time")
  par(new = TRUE)
  plot(1:n, smooth(sigma_female), type = "l", pch = 19, ylim = c(1, 3), xlim = c(1, n), col = "red", ylab = expression(sigma), xlab = "Time")
  legend("topleft", legend = c("Male", "Female"), col = c("Blue", "Red"), lty = c(1, 1), cex = 1.5)
  
  plot(1:n, smooth(gamma_male), type = "l", pch = 19, ylim = c(-0.2, 0), xlim = c(1, n), col = "blue", ylab = expression(xi), xlab = "Time")
  par(new = TRUE)
  plot(1:n, smooth(gamma_female), type = "l", pch = 19, ylim = c(-0.2, 0), xlim = c(1, n), col = "red", ylab = expression(xi), xlab = "Time")
  legend("topleft", legend = c("Male", "Female"), col = c("Blue", "Red"), lty = c(1, 1), cex = 1.5)
  
  #interaction effects 
  
}

LLGPD_repara_total <- function(paras, option){
  
  start_dates <- seq(as.Date("1811-01-01"), as.Date("1891-01-01"), "20 year")
  end_dates <- seq(as.Date("1831-01-01"), as.Date("1911-01-01"), "20 year")
  n <- length(start_dates)
  
  LL <- 0
  
  for (i in 1:n){
    
    start <- start_dates[i]
    end <- end_dates[i]
    grouped_data <- process_data(start, end,  u = 100)
    LL <- LL + LLGPD_repara(paras, grouped_data, option, u = 100, time = i)
    
    
    
  }
  
  return (-LL)
  

  
}


  
  
  


process_data <- function(start_date, end_date, u){
  
  data <- subset(uk, birth_year >= start_date & birth_year < end_date & Age >= u)
  grouped_data <- data
  grouped_data$Year <- NULL
  grouped_data$birth_year <- NULL
  grouped_data <- grouped_data %>% group_by(Age) %>% summarise_each(funs(sum))
  
  return (grouped_data)
  
}

test_dis <- function(){
  
  total_LL_GPD <- -opt_total_GPD$value
  total_LL_EXP <- -opt_total_EXP$value
  LR <- -2 * (total_LL_EXP - total_LL_GPD)
  p_value <- 1 - pchisq(LR, df = 1)
  
  return (p_value)
  
  
}

test_gender <- function(){
  
  LL_GPD_male <- -opt_male_GPD$value
  LL_GPD_female <- -opt_female_GPD$value
  total_LL_GPD <- -opt_total_GPD$value
  
  LR <- -2 * (total_LL_GPD - LL_GPD_female - LL_GPD_male)
  
  p_value <- 1 - pchisq(LR, df = 2)
  
  return (p_value)
  
}

LLGPD_repara <- function(paras, data, option, u, time){
  
  LL <- 0
  ages <- data$Age
  if (option == "male"){
    death <- data$Male
  }else if (option == "female"){
    death <- data$Female
  }else{
    death <- data$Total
  }
  n <- length(ages)
  
  alpha_1 <- paras[1]
  alpha_2 <- paras[2]
  xi <- alpha_1 + alpha_2 * time
  
  beta_1 <- paras[3]
  beta_2 <- paras[4]
  sigma <- beta_1 + beta_2 * time
  new_para <- c(sigma, xi)
  for (i in 1 : (n-1)){
    Age <- ages[i]
    LL <- LL + death[i] * log(GPD_cdf(new_para, Age + 1 - u) - GPD_cdf(new_para, Age - u))
  }
  LL <- LL + death[n] * log(1 - GPD_cdf(new_para, 110 - u))
  
  return (LL)
  
}



LL_GPD <- function(paras, data, option, u){
  
  LL <- 0
  ages <- data$Age
  if (option == "male"){
    death <- data$Male
  }else if (option == "female"){
    death <- data$Female
  }else{
    death <- data$Total
  }
  n <- length(ages)
  for (i in 1 : (n-1)){
    Age <- ages[i]
    LL <- LL + death[i] * log(GPD_cdf(paras, Age + 1 - u) - GPD_cdf(paras, Age - u))
  }
  LL <- LL + death[n] * log(1 - GPD_cdf(paras, 110 - u))
  
  return (-LL)
}

LL_EXP <- function(paras, data, option, u){
  
  LL <- 0
  ages <- data$Age
  if (option == "male"){
    death <- data$Male
  }else if (option == "female"){
    death <- data$Female
  }else{
    death <- data$Total
  }
  n <- length(ages)
  for (i in 1 : (n-1)){
    Age <- ages[i]
    LL <- LL + death[i] * log(EXP_cdf(paras, Age + 1 - u) - EXP_cdf(paras, Age - u))
  }
  LL <- LL + death[n] * log(1 - EXP_cdf(paras, 110 - u))
  
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