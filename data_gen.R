#Data and function set up

library(randomForest)
library(dplyr)

#Generate data
f.m<-function(x){2*51480*x^7*(1-x)^7}

meany <- function(x,z,u,beta1 = 1,f.m = function(x){2*51480*x^7*(1-x)^7},beta2 = 1){
  beta1*x+f.m(z)+beta2*u
  #f.m(z)+beta2*u
}

gen.data <- function(n,seed = 4,mean.y = meany,var.error=1){
set.seed(seed)
z <- runif(n,0,1)
x <- sapply(z,function(i) rnorm(1,(i-0.5)^2,1))
#Change from U ~ Norm(X,0.05^2) + Norm(X,0.05^2) to Norm(X,0.1^2) + Norm(X,0.1^2) add difference.
u <- runif(n,0,6)+sapply(x, function(i) rnorm(1,i,0.05^2))+sapply(z, function(i) rnorm(1,i,0.05^2))
#u <- runif(n,0,6)+sapply(x, function(i) rnorm(1,i,0.05^2))+sapply(z, function(i) rnorm(1,i,0.05^2))
y <- rnorm(n,mean.y(x,z,u),var.error)
data <- data.frame(y,x,u,z)
return(data)
}

#Generate missing
f.pis <- function(u,gamma0 = -2, gamma1 = 1, a1 = 0.5, a2 = 5.5){#missing model
  sta <- a1
  end <- a2
  temp <- exp(gamma0+gamma1*(u-sta)*(u>sta)*(u<=end)+gamma1*(end-sta)*(u>end))
  temp/(1+temp)
}

mis.data <- function(pis = f.pis,data = data){
  pi <- pis(data[,"u"])
  n = dim(data)[1]
  r <- rbinom(n,1,pi)
  data$r=r
  data$y[r==0] <- NA 
  return(data)
}
#f.m<-function(x){2*51480*x^7*(1-x)^7}
#t.data <- mis.data(data = gen.data(n=500))
#summary(3/5*(t.data$u-t.data$x-t.data$z-3)^2*(1/f.pis(t.data$u)-1))
# mean(1/f.pis(t.data$u)-1)
#model.mis <- mmis(t.data,name.outcome = "y")
