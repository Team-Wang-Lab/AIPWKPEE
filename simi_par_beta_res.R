  #Wrapper function for Estimating the final set of z.est for all z locations
  ###Note: Provide similar outcome of the beta estimates.
  semi.par.beta.res <- function(type ="aipw", x,y,z, h1a, h1b, zgrid = seq(min(z),max(z),length.out = 20), r,pis,delta,
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
    
    v.m <- est.var(y=y,x=x,z=z,f.m=z.init,data=data.frame(z=z,x=x),pi=pis,delta=delta,V="linear")#!!!!V can be changed to manual functions
    x.init <- betaest(zests = z.init,z=z[r==1], x= x[r==1],y=y[r==1],r=r[r==1],V.model=v.m,family=family,
                      delta=delta,type="naive",pis=pis,outcometype=outcometype)$estimate[1]
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
      v.m <- est.var(y=y,x=x,z=z,f.m=z.est,data=data,pi=pis,delta=delta,V="linear")
      beta.model <- betaest(zests = z.est, z=z,x= x,y=y,r=r,V.model=v.m,family=family,
                            delta=delta,type=type,pis=pis,outcometype=outcometype)
      beta.est <- beta.model$estimate[1]
      print(paste("iter=",iter))
      #Updata iteration time.
      iter = iter + 1
    }
    return(list(beta.est,z.est,beta.model,z.model))
  }
  
  
  
  #type ="aipw"; 
  #x=data$x;
  #y=data$y;z=data$z;
  #h1a=0.07;
  #h1b=0.15; 
  #zgrid = seq(min(z),max(z),length.out = 20); r=r
  #pis=pis;delta=delta;tol=1e-003;maxit=2;family=gaussian(link = "identity")