
GEV_LL <- function(paras, region){
  
  
  
  mu <- paras[1]
  sigma <- paras[2]
  xi <- paras[3]
  
  x <- region$total_age
  
  if (abs(xi) < 1e-9){
    
    LL <- -log(sigma) - exp(-(x - mu)/sigma) - (x - mu)/sigma
  }else{
    LL <- -log(sigma) - (1 + xi * (x - mu)/sigma) ^ (-1/xi) - (1 + 1/xi) * log(pmax(0, 1 + xi * (x - mu)/sigma))
    
  }
  
  
  
  return (LL)   
  
}




GEV_total_LL <- function(paras, region, k){
  
  LL <- 0
  region$total_age <- region$AGEYEARS + region$DAYSSINCEBD/365
  years <- as.character(seq(2000, 2017, 1))
  m <- length(years)
  for (i in 1:m){
    
    year <- years[i]
    region_year <- subset(region, year(DDATE) == year)
    region_year <- region_year[with(region_year, order(total_age)), ]
    n <- min(k, dim(region_year)[1])
    region_max <- tail(region_year, n)
    LL <- LL + GEV_LL_block(paras, region_max)
    
  }
  return (-LL)
  
}

GEV_LL_block <- function(paras, region){
  
  if (dim(region)[1] == 0){
    
    return (0)
  }
  
  mu <- paras[1]
  sigma <- paras[2]
  xi <- paras[3]
  
  region$total_age <- region$AGEYEARS + region$DAYSSINCEBD/365
  x <- region$total_age
  x <- sort(x)
  n <- length(x)
  inside <- pmax(0, (1 + xi * ((x -mu)/sigma)))

  
  LL <- -(inside[1])^(-1/xi) - n * log(sigma) - (1/xi + 1) * sum(log(inside))
  
  
  return (LL)
  
}
