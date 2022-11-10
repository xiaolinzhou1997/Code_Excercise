library(tidyverse)
library(hdm)
library(data.table)
library(furrr)
data(BLP)
BLP <- BLP$BLP
BLP$price <- BLP$price + 11.761
#Test dataset, cross-sectional.
BLP_test <- filter(BLP, cdid==10)[,7:11] %>%
  as.matrix()
BLP_share <- filter(BLP, cdid==10)[,"share"]
p <- filter(BLP, cdid==10)[,6]

{ #Parameters for testing
  sigma <- c(1,1,1,1,1)
  alpha <-  1
  delta <- runif(nrow(BLP_test),5,10)
}

n <- 1000
###Draw numbers for each year aka market.
log.y_draw <- randtoolbox::halton(n, 103, normal = TRUE,start = 20221108)*0.84+3.082
yp_draw <- sweep(exp(log.y_draw), 2, p) 
# log.yp <- matrix(rep(rep(log.y_draw), n),nrow=n)
# May want to use Nevo(2000) method to use y-p to capture price sensitivity. 
v_draw <- randtoolbox::halton(n, 5, normal = TRUE,start = 19970107)

v_draw <- split(v_draw, seq(nrow(v_draw)))
yp_draw <- split(yp_draw, seq(nrow(yp_draw)))
v_draw <- map(v_draw, as.matrix) %>% map(t)
yp_draw <- map(yp_draw, as.matrix)

cal_eu <- function(yp_draw, v_draw, char, alpha, delta, sigma){
  u <- delta + rowSums(char %*% sigma %*% v_draw) + alpha * t(yp_draw)
  eu <- u %>% exp()
    return(eu)
}



cal_share <- function(delta. = delta, ## J vector = mean utility
                     data,      ## J by K = from data where K = K characteristics
                     income,      ## S vector = income draw
                     yp,     ## S by J Matrix
                     v,      ## S by K = Var of K characteristics
                     p.=p,      ## J vector
                     alpha. = alpha,  ## scalar = coefficient of log.yp
                     sigma. = sigma)  ## K vector
{
  numer <- map2(yp, v, cal_eu, data, alpha., delta., sigma.)
  denom <- numer %>%
    map(sum) %>% 
    map(~ . + 1)
  share <- map2(numer, denom, ~ .x/.y) %>%
    map(as.data.table) %>%
    rbindlist() %>%
    colMeans()
  return(share)
}


### Note before we want can calculate share, we want to make sure log.yp is 
### defined, i.e. income is greater than price but not the case in the original
### BLP dataset. Solution: use Taylor expansion 

inner_loop <- function(delta, true_share, max_iter = 1000, tol = 10e-14){
  iter <- 0
  for (iter in iter:max_iter) {
    share <- cal_share(delta, BLP_test, yp=yp_draw, v=v_draw, p=p, alpha=alpha, sigma=sigma)
    #share <- share + runif(nrow(BLP_test),10e-5,10e-3)
    n_delta <- delta + log(true_share)- log(share)
    delta <- n_delta
    d <- mean(abs(true_share-share))
    if (d < tol){
      cat(paste("converged at", iter), fill=TRUE)
      break}
    else{
    iter <- iter + 1
    if(iter %% 10 == 0){
    cat(paste("iter", iter, "distance", d), fill=TRUE)
      }
    }
  }
  return(delta)
}


profvis(delta_new <- inner_loop(delta, BLP_share))


