### tidyverse for data cleaning/data.table/furrr for future mapping.
library(tidyverse)
library(data.table)
library(furrr)

### Get Nevo 2001 data from PyBLP and some housekeeping
{
  library(reticulate)
  use_python("/Library/Frameworks/Python.framework/Versions/3.11/bin/python3")
  pyblp <- import("pyblp")
  cereal <- fread(pyblp$data$NEVO_PRODUCTS_LOCATION)
  demographic <- fread(pyblp$data$NEVO_AGENTS_LOCATION)
  # add constant
  cereal <- cereal[,constant:=1] %>%
    relocate(constant, .before=prices)
  # add outside shares
  cereal <- cereal[,out_shares:=1-sum(shares),market_ids] %>%
    relocate(out_shares, .after=shares)
}

demog <- demographic[,c("income", "income_squared", "age", "child")] %>% 
  as.matrix()
X1 <- cereal[,c("constant", "prices", "sugar", "mushy")] %>%
  as.matrix()
X2 <- cereal[,c("constant", "prices", "sugar", "mushy")] %>%
  as.matrix()


k <- 4
N <- 20 #Draws for each market
J <- cereal[,list("J"=length(product_ids)),market_ids]

T <- nrow(J)

{ #Parameters for testing
  sigma <- c(0.377,1.848,0.004,0.081)
  alpha <-  -4
  delta <- log(cereal[,shares])-log(cereal[,out_shares])
}

v_draw <- randtoolbox::halton(n*T, k, normal = TRUE,start = 19970107)
pi <- c(3.09, 0, 1.186, 0, 16.6, -0.66, 0, 11.6,
        -0.193, 0, 0.03, 0, 1.468, 0, -1.5, 0)

cal_eu <- function(X2, DEMO, V, delta, sigma, pi, T, N){
  u <- matrix(nrow = T*24,ncol = N)
 for (t in 1:T) {
   for (n in 1:N) {
     for (j in 1:24) {
       u[(t-1)*24+j,n] <- delta[(t-1)*24+j] + sum(X2[(t-1)*24+j,] * V[(t-1)*20 + n,] * sigma) + 
         sum(rep(X2[(t-1)*24+j,],4) * rep(demog[(t-1)*20 + n,]) * pi)
     }
   }
 }
  eu <- u %>% exp()
    return(eu)
}


share <- 1e10

cal_share <- function(X2, demog, v_draw, delta, sigma, pi, T, N)  ## K vector
{
  utility <- cal_eu(X2, demog, v_draw, delta, sigma, pi, T, N)
  for (t in 1:T) {
    numer <- utility[(1+24*(t-1)):(24*t),] 
    denom <- 1 + colSums(utility[(1+24*(t-1)):(24*t),])
    share[(1+24*(t-1)):(24*t)] <- (t(numer)/denom) %>% colMeans()
  }
  return(share)
}


### Note before we want can calculate share, we want to make sure log.yp is 
### defined, i.e. income is greater than price but not the case in the original
### BLP dataset. Solution: use Taylor expansion 

inner_loop <- function(delta, true_share, max_iter = 1000, tol = 10e-14){
  iter <- 0
  for (iter in iter:max_iter) {
    share <- cal_share(X2, demog, v_draw, delta, sigma, pi, T, N)
    #share <- share + runif(nrow(BLP_test),10e-5,10e-3)
    n_delta <- delta + log(true_share)- log(share)
    delta <- n_delta
    d <- max(abs(true_share-share))
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


profvis(delta_new <- inner_loop(delta, cereal$shares))


