library(tidyverse)
library(hdm)
data(BLP)
BLP <- BLP$BLP
BLP$price <- BLP$price + 11.761
#Test dataset, cross-sectional.
BLP_test <- filter(BLP, cdid==10)[,7:11] %>%
  as.matrix()
BLP_share <- filter(BLP, cdid==10)[,"share"]
p <- filter(BLP, cdid==10)[,6]
x <- BLP_test
{ #Parameters for testing
  sigma <- seq(0,1,1/4)
  alpha <- - 0.5
  delta <- runif(nrow(BLP_test),5,10)
}
v <- v_draw
y <- y_draw
n <- 1000
###Draw numbers for each year aka market.
log.y_draw <- randtoolbox::halton(n, nrow(BLP_test), normal = TRUE,start = 20221108)*0.84+3.082
v_draw <- randtoolbox::halton(n, normal = TRUE,start = 19970107)
yp_draw <- sweep(exp(log.y_draw),2,p)
# log.yp <- matrix(rep(rep(log.y_draw), n),nrow=n)
# May want to use Nevo(2000) method to use y-p to capture price sensitivity. 
v_draw <- randtoolbox::halton(n, normal = TRUE,start = 19970107)

cal_share <- function(delta. = delta, ## J vector = mean utility
                     data,      ## J by K = from data where K = K characteristics
                     income,      ## S vector = income draw
                     yp,     ## S by J Matrix
                     v,      ## S by K = Var of K characteristics
                     p.=p,      ## J vector
                     alpha. = alpha,  ## scalar = coefficient of log.yp
                     sigma. = sigma)  ## K vector
{
  y <- y - p
  numer <- exp(delta + alpha + x %*% sigma %*% t(v) + alpha * t(y))
  denom <- 1 + rowSums(exp(delta + x %*% sigma %*% t(v) + alpha * t(y)))
  share <- mean(numer/denom)
  return(share)
}


### Note before we want can calculate share, we want to make sure log.yp is 
### defined, i.e. income is greater than price but not the case in the original
### BLP dataset. Solution: use Taylor expansion 

inner_loop <- function(delta, share, true_share, max_iter = 1000, tol = 10e-8){
  iter <- 1
  share <- runif(nrow(BLP_test),0.5,1)
  for (iter in iter:max_iter) {
    #share <- cal_share(delta,data,income,yp,v,p,alpha,sigma)
    share <- share + runif(nrow(BLP_test),10e-5,10e-3)
    n_delta <- delta + log(true_share)- log(share)
    delta <- n_delta
    d <- mean(abs(true_share-share))
    if (d > 3){
      print(paste("converged at", iter))
      break}
    else{
    iter <- iter + 1
    print(paste("iter", iter, "distance", d))}
  }
  return(delta)
}




