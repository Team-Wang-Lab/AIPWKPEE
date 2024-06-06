#semipar_est with fixed bandwith not changing over iterations.
#Wrapper function for Estimating the final set of z.est for all z locations

semi.par.est <- function(type ="aipw", x,y,z, bands, zgrid = seq(min(z),max(z),length.out = 20), r,pis,delta,
                         tol=1e-02,maxit=3,family=gaussian(link = "identity"),p=1){
  #Bands is a list of 2 elements 1=bands for the whole dataset 2=bands for the zgrid
  if(length(bands[[1]])!=length(z)&length(bands[[1]])!=1) stop("Length does not match for bandwiths and z.")
  if(length(bands[[2]])!=length(zgrid)&length(bands[[2]])!=1) stop("Length does not match for bandwiths and zgrid.")
  if(family$family=="gaussian"){
    m1 <- function(x) {x}
    outcometype="continuous"
  }else if(family$family=="binomial"){
    m1 <- function(x) {exp(x)/(1+exp(x))}
    outcometype="binary"
  }
  #Est starting beta(TRUE) and z.est(R==1) values. Using normal complete data without weights(Separately Now).
  z.init <- lpoly(type="naive", x=z[r==1],y=y[r==1],xpar = x[r==1],betapar = 2,h=0.2, xgrid = z[r==1], 
                 r=r[r==1],pis=pis[r==1],delta=delta[r==1],p=1)$f
  v.m <- est.var(y=y[r==1],x=x[r==1],z=z[r==1],f.m=z.init,data=data.frame(z=z[r==1],x=x[r==1]),pi=pis[r==1],delta=delta[r==1],V="linear",type=type)
  x.init <- betaest(zests = z.init,z=z[r==1], x= x[r==1],y=y[r==1],r=r[r==1],V.model=v.m,family=family,
                   delta=delta[r==1],type="naive",pis=pis[r==1],betastart=0,outcometype=outcometype)[1]
  print("Initiation done!")
  #Iterate Process
  beta.est <- x.init
  if(type=="aipw"){
  z.est <- rep(1,length(y))
  z.est[r==1] <- z.init
  z.old=rep(999,length(y))
  }else{
  z.est=z.init
  z.old=rep(999,length(y[r==1]))
  y=y[r==1]
  z=z[r==1]
  x=x[r==1]
  bands[[1]]<-bands[[1]][r==1]
  pis=pis[r==1]
  delta=delta[r==1]
  r=r[r==1]
  }
  beta.old=0
  iter=1
  while(iter <= maxit & mean(abs(z.old-z.est)/abs(z.est)) > tol & abs(beta.old-beta.est)/abs(beta.est) > tol){
    #Compare if the current iteration and the last iteration differs by what pct(Euclidian Now).
    print(max(abs(z.old-z.est)/abs(z.est)))
    print(abs(beta.old-beta.est)/abs(beta.est))
    #Est the z.est for all z locations.
    beta.old = beta.est 
    z.old = z.est 
    band.obj <- lpoly(x=z,y=y,xpar = x,betapar = beta.est,h=bands[[1]],xgrid=z,p=p,type=type,r=r,pis=pis,delta = delta)  
    z.est <- band.obj$f
    #Est beta corresponding to that.
    v.m <- est.var(y=y,x=x,z=z,beta=2,f.m=z.est,data=data,pi=pis,delta=delta,V="linear",type=type) 
    #v.m <- rep(1,length(y))
    #not close to the truth where intercept =0.91 but x is significant.
    beta.model <- betaest(zests = z.est, z=z,x= x,y=y,r=r,V.model=v.m,family=family,
                          delta=delta,type=type,pis=pis,betastart=beta.old,outcometype=outcometype)
    beta.est <- beta.model[1]
    beta.var <- beta.model[2]
    print(iter)
    #Updata iteration time.
    iter = iter + 1
    z.model.final <- lpoly(x=z,y=y,xpar = x,betapar = beta.est,h=bands[[2]],xgrid=zgrid,p=p,type=type,r=r,pis=pis,delta = delta)  
    z.est.final <- z.model.final$f
    z.est.var <-  z.model.final$var
  }
  return(list(beta.est=beta.est,z.est=z.est,beta.var=beta.var,z.est.final=z.est.final,z.est.var=z.est.var))
}

#########DATA EXPERIMENTS############
# mu <- function(t){#link function
#   t
# }
# 
# mu1 <- function(t){#Variance structure estimation.
#   rep(1, length(t))
# }
# 
# f.m<-function(x){2*51480*x^7*(1-x)^7}
# 
# f.theta <- function(z,beta0=1,beta1=1,mu.u=3){
#   beta0*f.m(z)+z+beta1*mu.u
#   #mu.u is the mean of U = (0+6)/2
# }
# 
# ns=500;seeds=5
# t.data <- mis.data(data = gen.data(n=ns,seed = seeds))
# data=t.data
# data = data %>% arrange(z)
# x <- data$x
# y <- data$y
# z <- data$z
# r <- data$r
# 
# p=1
# 
# f.theta.true <- function(xgrid) f.m(xgrid)+xgrid+3
# 
# delta <- f.theta(data$z)+2*data$x #Delta is true
# 
# pis = f.pis(data$u)#Pi model is true.
# 
# type ="aipw"
# h1a=0.1
# h1b=0.3
# zgrid = seq(0.05,0.95,length.out = 20)
# tol=5e-02
# maxit=3
# family=gaussian(link = "identity")
# p=1
# bands = list(rep(0.2,500),rep(0.2,length(zgrid)))
# ipw <- semi.par.est(x=x,z=z,y=y,type="ipw",bands = bands,zgrid = seq(0.05,0.95,length.out = 20),
#              tol=0.002,maxit = 5,p=1,r=r,pis=pis,delta=delta,family=gaussian(link = "identity"))
# 
# aipw <- semi.par.est(x=x,z=z,y=y,type="aipw",bands = bands,zgrid = seq(0.05,0.95,length.out = 20),
#              tol=0.002,maxit = 5,p=1,r=r,pis=pis,delta=delta,family=gaussian(link = "identity"))
# 
# mean(ipw$z.est.var-aipw$z.est.var)
