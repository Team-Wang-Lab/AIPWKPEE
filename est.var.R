# Estimate Variance

# Input: current f.m and beta; missing model; delta model; V structure; 
# Output: list(model, model coefficients, epsilon.star^2)
#First estimate f.m star
#Then model f.m star^2 with V structure(LINEAR OR A SPECIFIC FORMULA).

est.var <- function(y,x,z,data,f.m=function(x){2*51480*x^7*(1-x)^7},beta = 2,mis.model=NULL,pi,delta,V="linear",type="aipw"){
  #V(x,z)
  if(!(V=="linear"|class(V)=="formula")) stop("Wrong format for V. It should be eighter \"linear\" or a formula in character form.")
  r <- 1-is.na(y)
  y[is.na(y)]=0
  xfull=x
  zfull=z
  if(!is.null(mis.model)){
  pi <- predict(mis.model,newdata = data,type = "prob")[,2]#IS THAT PI FOR BEING OBSERVED OR MISSING?
  pi <- ifelse(pi==0,pi+0.001,pi)}
  if(class(f.m)=="function"){
  mu.hat <- f.m(z)+z+3+x*beta 
  } else {mu.hat <- f.m+x*beta } 
  if(type=="aipw"){
  epsilon.star <- as.vector(r/pi*(y-mu.hat)-(r/pi-1)*(delta-mu.hat ))
  }else if(type=="ipw"){
    pi=pi[r==1]
    z=z[r==1]
    x=x[r==1]
    y=y[r==1]
    mu.hat=mu.hat[r==1]
    r=r[r==1]
    epsilon.star <- as.vector(1/pi*(y-mu.hat))  
  }else{
    epsilon.star <- as.vector(y[r==1]-mu.hat[r==1])
    pi=pi[r==1]
    z=z[r==1]
    x=x[r==1]
  }
  epsilon.star=ifelse(epsilon.star==0,1,epsilon.star)
  if(V=="linear"){
    var.m <- lm(log(epsilon.star^2)~x+z)
    m.coef <- coef(var.m)
    names(m.coef) <- c("(Intercept)","x","z")
  }else{
    var.m <- lm(as.formula(paste("log(epsilon.star^2) ~",V,collapse = "")))
    m.coef <- coef(var.m)
  }
  fit.var = exp(predict(var.m,newdata=data.frame(x=xfull,z=zfull)))
  #Note the returning model is for log(var(Y|X,Z))
  return(list(var.model=var.m,model.coef=m.coef,fit.var = fit.var))
}

#Example
#data=t.data
#delta <- f.m(data$z)+2*data$x
#delta <- f.m(data$z)+3
#v.model <- est.var(y=data$y,x=data$x,z=data$z,data=data,pi=pis,delta=delta,V="linear")

#plot(v.model[[3]],10.232+2.189*data$x+4.876*data$z)
#v.model[[3]]
#fit.var <- 10.232+2.189*data$x+4.876*data$z
#cor(fit.var,v.model[[3]])#Correlation is only 0.19.
