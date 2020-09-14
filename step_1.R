setwd("/Users/august/Desktop/Dissertation")
complete_data <- read.csv("idl_complete_clean.csv")
data <- subset(complete_data, VALIDATION %in% c("Yes", "YES"))
data$DDATE<-as.Date(data$DDATE, "%d/%m/%Y")
data$BDATE<-as.Date(data$BDATE, "%d/%m/%Y")
data <- data[, 2:11]
data$BREGION <- NULL
data$DREGION <- NULL
threshold <- 105
data$achieve_threshold<- data$BDATE %m+% years(threshold)
data$total_age <- data$AGEYEARS + data$DAYSSINCEBD/365


set_threshold <- function(threshold){
  
  world <- subset(data, total_age >= threshold)
  world$excess_life <- world$total_age - threshold
  return (world)
  
}

#Deaths before 2017
world <- set_threshold(threshold)
N_Europe <- subset(world, DCOUNTRY %in% c("DNK", "DEU", "EW", "BEL", "FIN", "NOR", "SWE", "CHE"))
S_Europe <- subset(world, DCOUNTRY %in% c("ITA", "FRA", "ESP", "AUT"))
NN <- subset(N_Europe, DCOUNTRY != "EW")
SS <- subset(S_Europe, DCOUNTRY != "FRA")
FRA <- subset(world, DCOUNTRY == "FRA")
UK <- subset(world, DCOUNTRY == "EW")
CAN <- subset(world, DCOUNTRY == "CAN")
Europe <- rbind(N_Europe, S_Europe)


mortality <- function(paras, age){
  
  sigma <- paras[1]
  gamma <- paras[2]
  num <- ((1 + gamma * age/sigma)^(-1))
  
  return (num/sigma)
  
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

  return (-LL)
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

LLGPD_Unbiased_est <- function(region){
  
  est <- optim(c(1, 0.1), LLGPD_Unbiased, region = region)

  return (est$par)

}

LLGPD_Unbiased_est_CI <- function(region){
  
  est <- optim(c(1, 0.1), LLGPD_Unbiased, region = region)
  para <- est$par
  gamma <- para[2]
  hess <- hessian(LLGPD_Unbiased, para, region = region)
  sd.gamma <- sqrt(hess[2,2]/det(hess))
  CI <- c(gamma - 1.96*sd.gamma, gamma + 1.96*sd.gamma)
  
  return (CI)
}

LLEXP_Unbiased_est <- function(region){
  
  est <- optim(1, LLEXP_Unbiased, region = region, lower = 0.1, upper = 10, method = "Brent")
  
  return (est$par)
  
}

p_value_dis <- function(region){
  
  opt0 <- optim(c(1, 0.1), LLGPD_Unbiased, region = region)
  opt1 <- optim(1, LLEXP_Unbiased, region = region, lower = 0.1, upper = 10, method = "Brent")
  LL0 <- -opt0$value
  LL1 <- -opt1$value
  LR <- -2 * (LL1 - LL0)
  p <- 1 - pchisq(LR, df = 1)
  
  return (p)
  
}

p_val_gender <- function(region){
  
  female <- subset(region, SEX == "F")
  male <- subset(region, SEX == "M")
  
  opt_male_EXP <- optim(1, LLEXP_Unbiased, region = male, lower = 0.1, upper = 10, method = "Brent")
  opt_female_EXP <-  optim(1, LLEXP_Unbiased, region = female, lower = 0.1, upper = 10, method = "Brent")
  opt_total_EXP <-  optim(1, LLEXP_Unbiased, region = region, lower = 0.1, upper = 10, method = "Brent")
  LL_male_EXP <- -opt_male_EXP$value
  LL_female_EXP <- -opt_female_EXP$value
  LL_total_EXP <- -opt_total_EXP$value
  LR_EXP <- -2 * (LL_total_EXP - LL_female_EXP - LL_male_EXP)
  p_EXP <- 1 - pchisq(LR_EXP, df = 1)
  
  opt_male_GPD <- optim(c(1, 0.1), LLGPD_Unbiased, region = male)
  opt_female_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = female)
  opt_total_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = region)
  LL_male_GPD <- -opt_male_GPD$value
  LL_female_GPD <- -opt_female_GPD$value
  LL_total_GPD <- -opt_total_GPD$value
  LR_GPD <- -2 * (LL_total_GPD - LL_female_GPD - LL_male_GPD)
  p_GPD <- 1 - pchisq(LR_GPD, df = 2)

  p <- c(p_GPD, p_EXP)
  
  return (p)
  
  
}

P_val_time <- function(region){
  
  region <- region[with(region, order(BDATE)), ]
  time <- region$BDATE
  n <- length(time)
  mid <- time[ceiling(n/2) + 1]
  first_half <- subset(region, BDATE < mid)
  last_half <- subset(region, BDATE >= mid)
  
  opt_first_EXP <- optim(1, LLEXP_Unbiased, region = first_half, lower = 0.1, upper = 10, method = "Brent")
  opt_last_EXP <-  optim(1, LLEXP_Unbiased, region = last_half, lower = 0.1, upper = 10, method = "Brent")
  opt_total_EXP <-  optim(1, LLEXP_Unbiased, region = region, lower = 0.1, upper = 10, method = "Brent")
  LL_first_EXP <- -opt_first_EXP$value
  LL_last_EXP <- -opt_last_EXP$value
  LL_total_EXP <- -opt_total_EXP$value
  LR_EXP <- -2 * (LL_total_EXP - LL_first_EXP - LL_last_EXP)
  p_EXP <- 1 - pchisq(LR_EXP, df = 1)
  
  opt_first_GPD <- optim(c(1, 0.1), LLGPD_Unbiased, region = first_half)
  opt_last_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = last_half)
  opt_total_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = region)
  LL_first_GPD <- -opt_first_GPD$value
  LL_last_GPD <- -opt_last_GPD$value
  LL_total_GPD <- -opt_total_GPD$value
  LR_GPD <- -2 * (LL_total_GPD - LL_first_GPD - LL_last_GPD)
  p_GPD <- 1 - pchisq(LR_GPD, df = 2)
  
  return (c(p_GPD, p_EXP))
  
}

plot_4.2 <- function(region){
  
  female <- subset(region, SEX == "F")
  male <- subset(region, SEX == "M")
  opt_male_GPD <- optim(c(1, 0.1), LLGPD_Unbiased, region = male)
  opt_female_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = female)
  age <- seq(0, 15, 0.01)
  
  par_male <- opt_male_GPD$par
  par_female <- opt_female_GPD$par
  mortality_male <- mortality(par_male, age)
  mortality_male[mortality_male <= 0] <- Inf
  mortality_female <- mortality(par_female, age)
  mortality_female[mortality_female <= 0] <- Inf
  
  
  region <- region[with(region, order(BDATE)), ]
  time <- region$BDATE
  n <- length(time)
  mid <- time[ceiling(n/2) + 1]
  first_half <- subset(region, BDATE < mid)
  last_half <- subset(region, BDATE >= mid)
  
  opt_first_GPD <- optim(c(1, 0.1), LLGPD_Unbiased, region = first_half)
  opt_last_GPD <-  optim(c(1, 0.1), LLGPD_Unbiased, region = last_half)
  
  par_first <- opt_first_GPD$par
  par_last <- opt_last_GPD$par
  mortality_first <- mortality(par_first, age)
  mortality_first[mortality_first <= 0] <- Inf
  mortality_last <- mortality(par_last, age)
  mortality_last[mortality_last <= 0] <- Inf
  
  par(mfrow=c(1,2))
  plot(age + 105, mortality_male, type = "l",  ylim = c(0, 1), xlim = c(105, 120), col = "blue", ylab = "Force of Mortality", xlab = "Age")
  par(new = TRUE)
  plot(age + 105, mortality_female, type = "l", ylim = c(0, 100), xlim = c(105, 120), col = "red", ylab = "Force of Mortality", xlab = "Age", axes=FALSE)
  axis(4, ylim=c(0,100), col="red",col.axis="red",las=1)
  legend("topleft", col = c("blue", "red"), legend = c("Male (Left Axis)", "Female (Right Axis)"), lty = c(1, 1), cex = 1.3)
  
  plot(age + 105, mortality_first, type = "l", ylim = c(0, 2), xlim = c(105, 115), col = "blue", ylab = "Force of Mortality", xlab = "Age")
  par(new = TRUE)
  plot(age + 105, mortality_last, type = "l", ylim = c(0, 2), xlim = c(105, 115), col = "red", ylab = "Force of Mortality", xlab = "Age")
  legend("topleft", col = c("blue", "red"), legend = c("1856-1901", "1902-1909"), lty = c(1, 1), cex = 1.3)
  
  
}

plot_GPD <- function(){
  
  x <- seq(0.01, 5, 0.01)
  y1 <- dgpd(x, scale = 1, shape = -0.8)
  y2 <- dgpd(x, scale = 1, shape = -0.1)
  y3 <- dgpd(x, scale = 1, shape = 0.1)
  y4 <- dgpd(x, scale = 1, shape = 0.8)
  
  plot(x, y1, xlim = c(0, 5), ylim = c(0, 1), xlab = "X", ylab = "Probability Desnity", col = "skyblue", lty = 1, "l")
  par(new = TRUE)
  plot(x, y2, xlim = c(0, 5), ylim = c(0, 1), xlab = "X", ylab = "Probability Desnity", col = "deeppink", lty = 1, "l")
  par(new = TRUE)
  plot(x, y3, xlim = c(0, 5), ylim = c(0, 1), xlab = "X", ylab = "Probability Desnity", col = "orange", lty = 1, "l")
  par(new = TRUE)  
  plot(x, y4, xlim = c(0, 5), ylim = c(0, 1), xlab = "X",ylab = "Probability Desnity", col = "green", lty = 1, "l")
  legend("topright", col = c("skyblue","deeppink", "orange", "green"), legend = c(expression(paste(xi, " = ", -0.8)), expression(paste(xi, " = ", -0.1)),expression(paste(xi, " = ", 0.1)), expression(paste(xi, " = ", 0.8))), lty = c(1, 1, 1, 1), cex = 1.3)
  
  
}

plot_hist <- function(){
  
  par(mfrow=c(1,2))
  hist(UK$excess_life, breaks=seq(0, 20, 1), xlab = "Excess Age", ylab = "Count")
  hist(FRA$excess_life, breaks = seq(0, 20, 1), xlab = "Excess Age", ylab = "Count")
  
  
}