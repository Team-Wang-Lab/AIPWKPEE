#Functions for calculating pupurses

#mu.prime = return a vector of mu.prime based on outcome type.

mu.prime <- function(zest,beta,x,outcometype = "continuous"){#zest is theta(zi,beta) one single value. or a vector.
  if(length(zest) != length(x)) stop("Length of zest does not agree with length of x.")
  if(outcometype == "continuous"){
    mu.p <-rep(1,length(zest))   
  }else{
    mu.p <-exp(zest+x*beta)/(1+exp(zest+x*beta))
  }
  return(mu.p)
}


# phi AIPW/IPW

phi.aipw <- function(zest_i,z_i,
                     V.model,x,beta,outcometype = "continuous"){#V.model used zest vector as input for f.m!!!
  #mu.prime.zi <- sapply(1:length(x), function(i) mu.prime(zest=zest_i,beta=beta,x=x[i],outcometype = outcometype))
  mu.prime.zi <- mu.prime(zest=zest_i,beta=beta,x=x,outcometype = outcometype)
  #v.zi <- exp(predict(V.model$var.model,newdata = data.frame(z=z_i,x=x))) #zest_i is an single value. x should be a single value.
  v.zi <- V.model
  #print(v.zi)
  return(-mu.prime.zi^2/v.zi*x/mu.prime.zi^2/v.zi)
}


# epsilon star

eps.star <- function(pis,y,r,x,beta,zest,delta,outcometype = "continuous",type){
  if(outcometype=="continuous"){mu <- function(x) x
  }else{ mu <- function(x) exp(x)/(1+exp(x))}
  y <- ifelse(is.na(y),0,y)
  if(type=="aipw"){eps <- 1/pis*(r*(y-mu(x*beta + zest))-(r-pis)*(delta-mu(x*beta + zest)))
  }else if(type=="ipw"){eps <- 1/pis*(r*(y-mu(x*beta + zest)))
  }else if(type=="naive"){eps <- r*(y-mu(x*beta + zest))}
  return(eps)
}


#Var beta AIPW
var.beta<- function(zest,beta,x,z,pis,y,r,outcometype="continuous",V.model,
                    delta,type){#V.model used zest vector as input for f.m!!!
  mu.prime.vec <- mu.prime(zest=zest, beta=beta, x=x,outcometype = outcometype)
  V <- exp(predict(V.model$var.model,newdata = data.frame(z=z,x=x)))#V is modeling log 
  phi.aipw.vec <- sapply(1:length(z), function(i)
                  phi.aipw(zest_i=zest[i],z_i=z[i],V.model=V[i],x=x[i],beta=beta,outcometype = outcometype))#This could be wrong
  
  A.beta <- mean(mu.prime.vec^2*1/V*(x+phi.aipw.vec)^2)
  B.beta <- mean((mu.prime.vec/V*(x+phi.aipw.vec)*eps.star(pis,y,r,x,beta,zest,delta,outcometype,type))^2)
  n <- ifelse(type %in% c("ipw","aipw"),length(y),length(y[r==1]))
  #print(n)
  #n=length(y)
  betavar <- 1/A.beta * B.beta /A.beta/n
  return(betavar)
}


