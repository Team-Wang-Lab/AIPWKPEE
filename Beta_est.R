#Beta.est
betaest <- function(zests,z,x,y,r,V.model,
                    family=gaussian(link = "identity"),delta,type="aipw",pis,betastart,outcometype="continuous"){
  if(length(zests)!=length(x)) stop("Dimensions of theta.z and x does not match.")
  y <- c(y); x <- c(x); zests <- c(zests)
  n <- length(x)
  if(class(V.model)=="list"){
  V <- exp(predict(V.model$var.model,newdata = data.frame(z=z,x=x)))#V dimension is wrong.
  } else{
    V=V.model
  }
  y <- ifelse(is.na(y),0,y)
  
  data <- data.frame(x=x,y=y,zests=zests,r=r,delta=delta,type=type,pis=pis,V=V)
  #Z = design matrix for the non-parametric vector.
  
  beta_est_fun <- function(data){
    Y <- data$y
    X <- data$x
    zests=data$zests
    r=data$r
    delta=data$delta
    type=unique(data$type)[1]
    pis=data$pis
    V=data$V
    function(theta){
      mu = X*theta+zests
      if (type == "naive"){
        c(r *1/V*(Y - mu)*X)
      }else if (type == "ipw"){
        c(r / pis*1/V *(Y - mu)*X)
      } else if (type == "aipw"){
        c(1/V * (r / pis * (Y - mu) - (r / pis - 1) * (delta - mu))*X)
      }
    }
  }
  results <- m_estimate(
    estFUN = beta_est_fun,
    data=data,
    root_control = setup_root_control(start = c(betastart))
  )
  #varbeta <- var.beta(zest=zests,beta=results@estimates,x=x,z=z,pis=pis,y=y,r=r,outcometype="continuous",
  #                    V.model=V.model,delta=delta,type=type)
  return(c(results@estimates,results@vcov))
}

# 
# t.data <- mis.data(data = gen.data(n=500,seed = 7))
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
# pis = f.pis(data$u)
# zests <- f.m(z)+z+3
# delta <- 2*data$x+zests
# v.m <- est.var(y=y,x=x,z=z,f.m=zests,data=data,pi=pis,delta=delta,V="linear")
# 
# a<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
# family=gaussian(link = "identity"),delta=delta,type="aipw",pis=pis,betastart=2,outcometype="continuous")
# 
# b<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
#            family=gaussian(link = "identity"),delta=delta,type="ipw",pis=pis,betastart=2,outcometype="continuous")
# 
# c<-betaest(zests=zests,z=z,x=x,y=y,r=r,V.model=v.m,
#             family=gaussian(link = "identity"),delta=delta,type="naive",pis=pis,betastart=2,outcometype="continuous")
# a
# b
# c
# 
# 
# a;b;c
# sqrt(a[2])
# sqrt(b[2])
# sqrt(c[2])

