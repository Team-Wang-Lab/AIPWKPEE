#Zest
#Full version
#zest and lpoly for semi-parametric

#zest estimate 
#x = centered predictor in linear term (x_i-xgrid)/h
#y = outcome variable
#betastart = current estimate of alpha_0 and alpha_1 which is the partial derivative of f.theta at current alpha_0(how do I know f.theta?)
#weights = K_h(x_i-x)
#delta = Ymean
#r = indicator for observed
#pis = missing model probability of being observed
#mu = link function
zest <- function(x,y,betastart,xpar,betapar,weights, delta, r, pis, type,mu=function(x){x}){
  y <- c(y); x <- c(x); betastart <- c(betastart); weights <- c(weights)
  y <- ifelse(is.na(y),0,y)
  p <- length(betastart) - 1 #p is the polynomial term power
  n <- length(x)
  X <- rep(1,n)
  for (i in 1:p) { X <- cbind(X, x^i) }#X design matrix = [1 x_i] for p = 1
  
  beta <- betastart
  
  #Z = design matrix for the non-parametric vector.
  
  z.fun <- function(X, y,xpar,betapar, weights, delta, r, pis, type,beta){
    #weights=Kernel func
    #delta=the delta model delta(X,Z,U)
    #r=1 - I(Y_i=missing)
    #pis = weights = prob of being observed
    if (type == "naive"){
      #print(X%*%beta)
      return(weights * r *(y - mu(X%*%beta+xpar*betapar))* X)
    }
    if (type == "ipw"){
      return(weights * r / pis *(y - mu(X%*%beta+xpar*betapar)) * X)
    }
    if (type == "aipw"){
      return(weights * (r / pis * (y - mu(X%*%beta+xpar*betapar)) - 
                          (r / pis - 1) * (delta - mu(X%*%beta+xpar*betapar))) * X)
    }
  }
  objective <- function(beta){
    Z.final <- rowSums(sapply(1:n, function(i) z.fun(X[i,], y[i],xpar[i],betapar ,weights[i], delta[i], r[i], pis[i], type,beta)))
    #creates many NAs -- change NA y's to 0
    Z.final %*% Z.final
  }
  
  beta <- optim(betastart, objective)[[1]]
  
  #Calculate variance
  calcvar <- function(X, y,xpar,betapar, delta, pis, type,beta){
    
    if (type == "ipw"){
      return(3/5*(1/pis+(delta-mu(X%*%beta+xpar*betapar))^2)/pis)
    }
    if (type == "aipw"){
      return(3/5*(1/pis+(delta-mu(X%*%beta+xpar*betapar))^2)) 
    }
  }
  
  if(type=="aipw"){
    #varz = mean(sapply(1:n,function(u) 
    #  3/5*(1/pis[u]+(delta[u]-mu(sum(X[u,]*beta)+xpar[u]*betapar))^2)))
    varz <- rowMeans(matrix(sapply(1:n, function(i) calcvar(X[i,], y[i],xpar[i],betapar , delta[i], pis[i], type,beta)),nrow = 1,byrow = TRUE))
    print(varz)
  }else if(type=="ipw"){
    #varz = mean(sapply(1:n,function(u) 3/5*(1+(delta[u]-mu(sum(X[u,]*beta)+xpar[u]*betapar))^2)/pis[u]))  
    varz <- rowMeans(matrix(sapply(1:n, function(i) calcvar(X[i,], y[i],xpar[i],betapar , delta[i], pis[i], type,beta)),nrow = 1,byrow = TRUE))
  }else{
    W <- mu1(X%*%beta)# W is the weight matrix all 1 s
    W.half <- sqrt(W)
    W <- diag(as.vector(W))
    W.half <- diag(as.vector(W.half))
    
    H <- (W.half %*% X) %*% solve(t(X) %*% W %*% X) %*% (t(X) %*% W.half)#H matrix asumming weight to be 1 across all.
    H <- diag(H)#h_{ii}?
    
    b.fun <- function(X, betapar,xpar,y, weights, delta, r, type){
      X <- matrix(X, ncol = 1)
      function(beta){
        return((weights * r * (y - mu(X%*%beta+betapar*xpar)))^2 * X %*% t(X))#weights = K_h(z_i-z)
      }
    }
    
    B <- matrix(0, nrow = 1+p, ncol = 1+p)
    for (i in 1:n){
      B <- B + b.fun(X[i, ], betapar,xpar[i],y[i], weights[i], delta[i], r[i], type)(beta) / (1 - H[i])#?????????????
    }
    B <- B/n
    
    a.fun <- function(X, betapar,xpar,y, weights, delta, r, type){
      X <- matrix(X, ncol = 1)
      function(beta){
        return(weights * r * mu1(X%*%beta+betapar*xpar) * X %*% t(X))
      }
      
    }
    A <- matrix(0, nrow = 1+p, ncol = 1+p)
    for (i in 1:n){
      A <- A + a.fun(X[i, ], betapar,xpar[i],y[i], weights[i], delta[i], r[i], type)(beta)
    }
    A <- A/n
    f <- function(m) class(try(solve(m),silent=T))=="matrix"
    
    if(f(A)==FALSE) stop("A does not have full rank.")
    Ainv <- solve(A)
    varz <- Ainv %*% B %*% Ainv
  }
  return(list(beta=beta, varmatrix=varz))
}


### Estimate theta(z_i,beta) on a grid and a given bandwidth. ####
# lpoly output theta(Z_i,beta) for xgrid or Z_i and the variances of the estimate.
# x = a vector of a variable that has non-linarly related to y(z in our manuscript).
# x1 = a vector of a variable that is linearly related to y
# y = a vector of outcome variable. Can have missingness.
# xgrid = the grid for variable x that we will evaluate theta(x_i,beta) on.
# p = order of the polynomial.
# type = c("naive","ipw","aipw")
# r = indicator of observed index in y. Can leave empty. 
# pis = a vector of probability for being observed.
# delta = estimated mean of y.
## Output: 
# f = theta(x_i,beta) = alpha_0
# var = var(alpha_0)


#lpoly is so slow now. what cause this?
lpoly <- function(x,y,xpar,betapar,h,xgrid,p,
                  type = type,r = NULL,pis = pis,delta = delta){
  # vectors:
  x<-c(x); y<-c(y); 
  if(is.null(r)) r = 1 - is.na(y)
  h<-c(h); xgrid<-c(xgrid)
  
  n <- length(x)
  m <- length(xgrid)
  
  if (length(h) == 1) h <- rep(h,m) 
  #single bandwith
  
  f <- rep(0,m) # length of xgrid.
  var <- f # 
  # initiate the betas p dimension.
  betastart = rep(1,p+1)
  #v.mz <- est.var(y=y,x=xpar,z=x,beta=betapar,f.m=function(x){2*51480*x^7*(1-x)^7},data=data.frame(z=x,x=xpar),pi=pis,delta=delta,V="linear")
  for (i in 1:m){#iteratively update across xgrid
    xch <- (x - xgrid[i]) / h[i]
    #This h[i] should be in the kernel only but now it is on the centered x. This might lead to the difference in zest.
    w <- (3/4)*(1 - (xch)^2)*(abs(xch) < 1)#this is the kernel function [plug in xc to kernel] K_h(z_i-z)
    xc <- (x - xgrid[i])
    
    out <- zest(x=xc,y=y,xpar=xpar,betapar=betapar,betastart=betastart,weights=w,type = type,r=r,pis=pis,delta = delta) 
    #Variance are negative.
    
    # Why betastart has two values? one is the f.theta(current xgrid) the other is partial f.theta at xgrid.
    # out - a list of two elements: beta estimate, variance covariance matrix of the beta estimate.
    beta <- out$beta
    varmatrix <- out$varmatrix / n #This is the variance of sigma(z_i) original variance is the variance of sqrt(n)(hat.sigma-sigma)
    
    f[i] <- beta[1]
    if(type=="naive"){
      var[i] <- varmatrix[1,1]
    }else{
      var[i]=varmatrix
      #print(var[i])
    }
    betastart <- beta
  }
  
  return(list(f=f,var=var))
}


ebbs <- function(fun, x,y,xpar,betapar,h1a,h1b,p=1,msespan=0,M1=14,J1=1,J2=2,nterms=2,
                 xgrid=seq(min(x), max(x), length=20),nskip=0,bandspan=4,
                 type,r,pis,delta,gridallx=TRUE){
  
  ##fix the grid for bandwith h.
  h1a=log(h1a)
  h1b=log(h1b) 
  m=length(xgrid)
  
  #If it is single bandwith for all z_i
  if (length(h1a) == 1) h1a = rep(h1a, m)
  if (length(h1b) == 1) h1b = rep(h1b, m)
  h2 <- matrix(sapply(1:m,function(i)exp(seq(h1a[i],h1b[i], length=M1))),nrow = m, ncol = M1,byrow = TRUE)
  
  #What is this?  Find nearby points to estimate bias at one bandwidth.  
  maxJ10 = max(0, J1)
  fvect = matrix(0, nrow = m, ncol = M1)#m x M1
  var2 = fvect
  band = rep(0, m)
  msevect = matrix(0, nrow = m, ncol = M1-J2-maxJ10) #m x (M1-J2-maxJ10)
  msehat = rep(0, m)
  
  #Estimate alpha0 = theta(x,beta) under different bandwith(h).
  result=sapply(1:M1,function(k){ #M1 = different h's to pick
    #fun is the input function nomissing, naive ipw, aipw
    if(fun=="nomissing"){
      result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=h2[, k], xgrid=xgrid, p=p, delta=NA, r=rep(1, length(y)), pis=NA, type="naive")
    }else if(fun=="naive"){
      result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=h2[, k], xgrid=xgrid, p=p, delta=NA, r=r, pis = NA,type = "naive")
    }else if(fun=="ipw"){
      result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=h2[, k], xgrid=xgrid, p=p, delta=delta, r=r, pis=pis,type = "ipw")
    }else if(fun=="aipw"){
      result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=h2[,k], xgrid=xgrid, p=p, delta=delta, r=r, pis=pis,type = "aipw")
    }
    #result=fun(x,y,h2[, k],xgrid,p,...) #output alpha0 = theta(x,beta) and its variance 
    c(result$f,result$var) #alpha0's#variances
  })
  halfrow <- dim(result)[1]/2
  #print(dim(result))
  fvect <- result[1:halfrow,]
  var2 <- result[(halfrow+1):dim(result)[1],]
  #What is H for?
  H = matrix(1, nrow = J2 + J1 + 1, ncol = nterms+1) #I_{(J2+J1) x (nterms+1)}
  
  for (l in c(1:m)){#??????????? calculate coefficient bias.
    for (ll in c(( 1 + maxJ10 ):(M1-J2))){
      bot = ll-J1
      top = ll+J2
      #print(c(bot,top))
      for (j in c(2:(1+nterms))){
        H[, j] = h2[l, bot:top] ^ (p+j-1)
      }
      biascoef = solve(t(H) %*% H + .0000001*diag(nterms+1)) %*% ( t(H) %*% fvect[l, (bot:top)] )
      h3= h2[l,ll] ^ ( (p+1):(p+nterms) )
      biassqu = ( biascoef[2:(nterms+1)] %*% h3 )^2  
      msevect[l,ll-maxJ10] = biassqu + var2[l,ll]
    }
  }
  
  n2=35
  msevect2=matrix(1, nrow = m, ncol = n2) # INTERPOLATE MSE TO FINER GRID
  
  eps <- .Machine$double.eps
  h4=matrix(0, nrow = m, ncol = n2)
  
  for (l in 1:m){
    h4[l,] <- seq(h2[l,1+maxJ10]+2*eps,h2[l,M1-J2]-2*eps,length=n2)
    msevect2[l, ] = interp1(h2[l,(maxJ10+1):(M1-J2)],msevect[l, ],h4[l, ], "cubic")#the points are all NA.
  }
  
  kopt = 2
  if (msespan == 0){
    mseweight = 1
  } else mseweight = c(1:(msespan+1), msespan:1) # WEIGHTS FOR SMOOTHING THE MSE
  
  mseweight = mseweight / sum(mseweight)
  
  for (l in (1+msespan):(m-msespan)){
    index = (l-msespan):(l+msespan)
    mse = matrix(mseweight, nrow = 1) %*% matrix(msevect2[index, ], nrow = 1 + 2*msespan)
    im=  n2
    for (lll in 1:(n2-kopt)){
      if (min(mse[(lll+1):(lll+kopt)]) > mse[lll]){
        if (im == n2){
          im =lll
        }
      }
    }
    band[l] = h4[l,im]
    msehat[l] = mse[im]
  }
  
  msns = msespan + nskip #  ESTIMATE MSE NEAR BOUNDARIES
  if (msns > 0){
    band[1:msns]   = rep(band[msns +1], msns)
    band[(m-msns+1):m] = rep(band[m-msns], msns)
  }
  
  
  if (bandspan > 0){ #   %  SMOOTH THE BANDWIDTH Creates NAN: yes if the min(xgrid) is less than the min(x) then it might create NA.
    band2=band 
    delta1 = bandspan * (xgrid[2] - xgrid[1])#delta is the distance of neighborhood = bandspan * diff in the first two elements
    for (i in 1:m){
      #If distance between the target and all others < delta we weight them based on distances.
      #If distance is greater than delta, then we think it is outside of the neighborhood.
      wt = pmax(delta1 - abs(xgrid-xgrid[i]),rep(0, m))
      wt = wt / sum(wt)
      band2[i] = wt %*% band
    }		
    band = band2
  }
  xgridold = xgrid 
  #TRY IF SET THIS INTO ACTUALL Z VS SETTING THIS TO THE 20 POINTS. COMPARE RESULTS.
  if(gridallx==TRUE){
    xgrid = x
    band = interp1(xgridold,band,xgrid,"linear")
  }
  band = mean(band)
  if(fun=="nomissing"){
    result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=band, xgrid=xgrid, p=p, delta=NA, r=rep(1, length(y)), pis=NA, type="naive")
  }else if(fun=="naive"){
    result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=band, xgrid=xgrid, p=p, delta=NA, r=r, pis = NA, type = "naive")
  }else if(fun=="ipw"){
    result = lpoly(x=x, y=y, xpar = xpar,betapar = betapar,h=band, xgrid=xgrid, p=p, delta=delta, r=r, pis, type = "ipw")
  }else if(fun=="aipw"){
    result = lpoly(x=x, y=y,xpar = xpar,betapar = betapar, h=band, xgrid=xgrid, p=p, delta=delta, r=r,pis=pis, type = "aipw")
  }
  list(f = result$f, var = result$var, xgrid = xgrid, band = band)
}

interp1 <- function(X, Y, XI, method){
  if (method == "cubic")
    return(spline(X, Y, xout = XI)$y)
  if (method == "linear")
    return(approx(X, Y, xout = XI)$y)
}

band.z <- function(ebbs.obj,z){
  bandind <- sapply(z, function(i) {
    temp <- which(order(c(i,ebbs.obj$xgrid))==1)
    ifelse(temp>length(z),temp-1,temp)
  })
  return(ebbs.obj$band[bandind])
}

