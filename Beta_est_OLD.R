#Beta.est
betaest <- function(zests,z,x,y,r,V.model,
                    family=gaussian(link = "identity"),delta,type="aipw",pis,betastart,outcometype="continuous"){
  if(length(zests)!=length(x)) stop("Dimensions of theta.z and x does not match.")
  y <- c(y); x <- c(x); zests <- c(zests)
  n <- length(x)
  V <- exp(predict(V.model$var.model,newdata = data.frame(z=z,x=x)))#V dimension is wrong.
  beta.fun <- function(x, y, delta, r, pis, type,beta){
    #delta=the delta model delta(X,Z,U)
    #r=1 - I(Y_i=missing)
    #pis = 1/weights
    y <- ifelse(is.na(y),0,y)
    if (type == "naive"){
      return( r * 1/V*(y - mu(x * beta + zests))*x)
    }
    if (type == "ipw"){
      return( r / pis *1/V*(y - mu(x * beta + zests))*x)
              #r / pis *(y - mu(x * beta + zests))*x)#add partial x parital theta
    }
    if (type == "aipw"){
      return( 1/V*(r / pis * (y - mu(x * beta+zests)) - (r / pis - 1) * (delta - mu(x * beta+zests)))*x)
        #(r / pis * (y - mu(x * beta+zests)) - (r / pis - 1) * (delta - mu(x * beta+zests)))*x)
    }
  }
  objective <- function(beta){
    Beta.final <- rowSums(sapply(1:n, function(i) 
      beta.fun(x=x[i], y=y[i],  delta=delta[i], r=r[i], pis=pis[i], type=type,beta=beta)))
    t(Beta.final) %*% Beta.final
    #creates many NAs -- change NA y's to 0
  }
  
  beta <- optimize(f=objective,lower = betastart-2, upper=betastart+2,tol = 0.0001)$minimum 
    #optim(betastart,fn=objective,method = "Brent",lower = -10000,upper = 10000)[[1]]  
  
  varbeta <- var.beta(zest=zests,beta=beta,x=x,z=z,pis=pis,y=y,r=r,outcometype="continuous",
                      V.model=V.model,delta=delta,type=type)
  return(c(beta,varbeta))
}


# t.data <- mis.data(data = gen.data(n=500,seed = 4))
# data=t.data
# x <- data$x
# y <- data$y
# z <- data$z
# r <- data$r
# ngrids <- 20
# xgrid <- seq(min(z),max(z),length.out = 20)
# h=0.07
# p=1
# x=data$x
# y=data$y
# r=data$r
# z=data$z
# delta <- 2*data$x
# pis = f.pis(data$u)
# zests <- f.m(z)+z+3
# v.m <- est.var(y=y,x=x,z=z,f.m=zests,data=data,pi=pis,delta=delta,V="linear")
# 
# a<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
# family=gaussian(link = "identity"),delta=delta,type="aipw",pis=pis,betastart=1,outcometype="continuous")
# 
# b<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
#            family=gaussian(link = "identity"),delta=delta,type="ipw",pis=pis,betastart=1,outcometype="continuous")
# 
# c<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
#             family=gaussian(link = "identity"),delta=delta,type="naive",pis=pis,betastart=1,outcometype="continuous")
# 
# 
# a;b;c
# sqrt(a[2])
# sqrt(b[2])
# sqrt(c[2])

