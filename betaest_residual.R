#Beta estimate residual

betaest <- function(zests, z,x,y,r,V.model,family=gaussian(link = "identity"),delta,type="aipw",pis,outcometype){
  if(length(zests)!=length(x)) stop("Dimensions of theta.z and x does not match.")
  yori=y
  xori=x
  zestsori = zests 
  if(type=="aipw"){
    y[r==0] <- delta[r==0]
    y[r==1] <- (y[r==1]-delta[r==1])/pis[r==1]+delta[r==1]
    weights <- rep(1,length(x))
  }else if(type %in% c("ipw")){
    y <- y[r==1]
    x <- x[r==1]
    zests <- zests[r==1]
    weights <- 1/pis[r==1]
  }else if(type=="naive"){
    y <- y[r==1]
    x <- x[r==1]
    zests <- zests[r==1]
    weights <- rep(1,length(x))
  }
    beta.model <- glm(y~x, offset = zests,weights = weights,family = family)
    beta <- beta.model$coefficients[2] 
    beta.var <- (summary(beta.model)$coefficients[2,2])^2
    print(beta.var)
    varbeta <- var.beta(zest=zestsori,beta=beta,x=xori,z=z,pis=pis,y=yori,r=r,outcometype=outcometype,
                        V.model=V.model,delta=delta,type=type)
    return(list(beta.model=beta.model,estimate = c(beta,varbeta)))
}

#1. variance calculation after estimation. Plug in to get true variance.
#2. stopping criteria max(theta(z)-hattheta(z)/theta(z)) or (beta-hatbeta)/beta deviance < 0.01.[DONE]
#3. write beta.est into optimization problem by optium. 
set.seed(4)
t.data <- mis.data(data = gen.data(n=500,seed = 1234))
data=t.data
x <- data$x
y <- data$y
z <- data$z
r <- data$r
ngrids <- 20
xgrid <- seq(min(z),max(z),length.out = 20)
h=0.07
p=1
x=data$x
y=data$y
r=data$r
z=data$z
delta <- 2*data$x
pis = f.pis(data$u)
zests <- f.m(z)+z+3
v.m <- est.var(y=y,x=x,z=z,f.m=zests,data=data,pi=pis,delta=delta,V="linear")

a<-betaest(zests = zests, z=z,x= x,y=y,r=r,V.model = v.m,family=gaussian(link = "identity"),
           delta=delta,type="naive",pis=pis,outcometype = "continuous")
betaest(zests = zests, x= x,y=y,r=r,family=gaussian(link = "identity"),delta=delta,type="ipw",pis=pis)
betaest(zests = zests, x= x,y=y,r=r,family=gaussian(link = "identity"),delta=delta,type="naive",pis=pis)



semi.par.est <- function(type ="aipw", x,y,z, h1a, h1b, zgrid = seq(min(z),max(z),length.out = 20), r,pis,delta,
                         tol=1e-03,maxit=3,family=gaussian(link = "identity"),p=1){
  if(family$family=="gaussian"){
    m1 <- function(x) {x}
    outcometype="continuous"
  }else if(family$family=="binomial"){
    m1 <- function(x) {exp(x)/(1+exp(x))}
    outcometype="binary"
  }
  #Est starting beta and z.est values. Using normal complete data without weights(Separately Now).
  z.init <- lpoly(type="naive", x=z[r==1],y=y[r==1],xpar = x[r==1],betapar = 0,h=(h1b+h1a)/2, xgrid = z[r==1], 
                  r=r,pis=pis,delta=delta,p=1)$f
  
  #v.m <- est.var(y=y,x=x,z=z,f.m=z.init,data=data.frame(z=z,x=x),pi=pis,delta=delta,V="linear")#!!!!V can be changed to manual functions
  x.init <- betaest(zests = z.init,x= x[r==1],y=y[r==1],r=r[r==1],
                    delta=delta,type="naive",pis=pis)$estimate[1]
  print("Initiation done!")
  #Iterate Process
  #Est the bandwidth
  beta.est <- x.init
  z.est <- rep(0,length(y))
  z.est[r==1] <- z.init
  beta.old=0
  z.old=rep(999,length(y))
  iter=1
  while(iter <= maxit & max(abs(z.old-z.est)/abs(z.est)) > tol & abs(beta.old-beta.est)/abs(beta.est) > tol){
    #dist(rbind(old=c(beta.old,z.old),cur=c(beta.est,z.est)),method = "euclidian") > tol){#Compare if the current iteration and the last iteration differs by what pct(Euclidian Now).
    print(max(abs(z.old-z.est)/abs(z.est)))
    print(abs(beta.old-beta.est)/abs(beta.est))
    #Est the z.est for all z locations.
    beta.old = beta.est 
    z.old = z.est 
    #band.obj <- ebbs(fun=type, x=z,y=y-x*beta.old,h1a=h1a, h1b=h1b, xgrid = zgrid, r=r,pis=pis,delta=delta)
    band.obj <- ebbs(fun=type, x=z,y=y,xpar = x,betapar = beta.old,h1a=h1a, h1b=h1b, xgrid = zgrid, r=r,pis=pis,delta=delta)
    bands <- band.z(band.obj,z=z)
    z.model <- lpoly(x=z,y=y,xpar = x,betapar = beta.old,h=bands,xgrid=z,p=p,type=type,r=r,pis=pis,delta = delta)
    z.est <- z.model$f
    #Est beta corresponding to that.
    #beta.model <- betaest(zests = z.est, x= x,y=y,r=r,family=family,
    #                    delta=delta,type="naive",pis=pis)
    #v.m <- est.var(y=y,x=x,z=z,f.m=z.est,data=data,pi=pis,delta=delta,V="linear")
    beta.model <- betaest(zests = z.est, x= x,y=y,r=r,
                          delta=delta,type=type,pis=pis)
    beta.est <- beta.model$estimate[1]
    print(iter)
    #Updata iteration time.
    iter = iter + 1
  }
  return(list(beta.est,z.est,beta.model,z.model))
}


