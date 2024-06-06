#Wrapper function for Estimating the final set of z.est for all z locations

semi.par.est <- function(type ="aipw", x,y,z, h1a, h1b, zgrid = seq(min(z),max(z),length.out = 20), r,pis,delta,
                         tol=5e-02,maxit=3,family=gaussian(link = "identity"),p=1){
  if(family$family=="gaussian"){
    m1 <- function(x) {x}
    outcometype="continuous"
  }else if(family$family=="binomial"){
    m1 <- function(x) {exp(x)/(1+exp(x))}
    outcometype="binary"
  }
  #Est starting beta(TRUE) and z.est(R==1) values. Using normal complete data without weights(Separately Now).
  #z.init <- lpoly(type="naive", x=z[r==1],y=y[r==1],xpar = x[r==1],betapar = 2,h=(h1b+h1a)/2, xgrid = z[r==1], 
  #               r=r,pis=pis[r==1],delta=delta[r==1],p=1)$f
  z.init <- 2*51480*z^7*(1-z)^7+z+3
  #v.m <- est.var(y=y[r==1],x=x[r==1],z=z[r==1],f.m=z.init,data=data.frame(z=z[r==1],x=x[r==1]),pi=pis[r==1],delta=delta[r==1],V="linear")
  #!!!!V can be changed to manual functions not the same length z.init = 300 ish
  
  #x.init <- betaest(zests = z.init,z=z[r==1], x= x[r==1],y=y[r==1],r=r[r==1],V.model=v.m,family=family,
  #                 delta=delta,type="naive",pis=pis,betastart=0,outcometype=outcometype)[1]
  x.init <- 2
  print("Initiation done!")
  #Iterate Process
  #Est the bandwidth
  beta.est <- x.init
  #z.est <- rep(0,length(y))
  #z.est[r==1] <- z.init
  z.est=z.init
  beta.old=0
  z.old=rep(999,length(y))
  iter=1
  while(iter <= maxit & mean(abs(z.old-z.est)/abs(z.est)) > tol & abs(beta.old-beta.est)/abs(beta.est) > tol){
  #dist(rbind(old=c(beta.old,z.old),cur=c(beta.est,z.est)),method = "euclidian") > tol){
    
    #Compare if the current iteration and the last iteration differs by what pct(Euclidian Now).
    print(max(abs(z.old-z.est)/abs(z.est)))
    print(abs(beta.old-beta.est)/abs(beta.est))
    #Est the z.est for all z locations.
    beta.old = beta.est 
    z.old = z.est 
    
    if(iter==1){
    band.obj <- ebbs(fun=type, x=z,y=y,xpar = x,betapar = 2,h1a=h1a, h1b=h1b, xgrid = zgrid, r=r,pis=pis,delta=delta,gridallx=TRUE)
    bands = band.obj$band
    }else{
    band.obj <- lpoly(x=z,y=y,xpar = x,betapar = 2,h=bands,xgrid=z,p=p,type=type,r=r,pis=pis,delta = delta)  
    }
    #bands <- band.z(band.obj,z=z)
    #z.model <- lpoly(x=z,y=y,xpar = x,betapar = beta.old,h=bands,xgrid=z,p=p,type=type,r=r,pis=pis,delta = delta)
    #z.est <- z.model$f
    z.est <- band.obj$f
    #Est beta corresponding to that.
    #VM and betaest seems to be the real problem---- betaest is not good.
    v.m <- est.var(y=y,x=x,z=z,beta=2,f.m=z.est,data=data,pi=pis,delta=delta,V="linear") 
    #not close to the truth where intercept =0.91 but x is significant.
    beta.model <- betaest(zests = z.est, z=z,x= x,y=y,r=r,V.model=v.m,family=family,
                        delta=delta,type=type,pis=pis,betastart=beta.old,outcometype=outcometype)
    beta.est <- beta.model[1]
    beta.var <- beta.model[2]
    print(iter)
    #Updata iteration time.
    iter = iter + 1
    #eval.band <- band.z(band.obj,z=seq(0.05,0.95,length.out = 20))
    z.model.final <- ebbs(fun=type, x=z,y=y,xpar = x,betapar = 2,h1a=h1a, h1b=h1b, xgrid = seq(0.05,0.95,length.out = 20), r=r,pis=pis,delta=delta,gridallx=FALSE)
    #  lpoly(x=z,y=y,xpar = x,betapar = beta.old,h=eval.band,xgrid=seq(0.05,0.95,length.out = 20),
    #                       p=p,type=type,r=r,pis=pis,delta = delta)
    z.est.final <- z.model.final$f
    z.est.var <-  z.model.final$var
  }
  print(z.model.final$band)
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
# ns=500;seeds=12
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
# zgrid = seq(min(z),max(z),length.out = 20) 
# tol=5e-02
# maxit=3
# family=gaussian(link = "identity")
# p=1
