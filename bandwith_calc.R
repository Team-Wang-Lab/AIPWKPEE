int<-integrate(function(x){(-1 + x)^5*x^5*(3 - 13*x + 13*x^2)}, 0,1)
log((-1441440)*(int$value))

exp((log(4)+log(0.3^2)-log(1/5)-log((-1441440)*(int$value)))/5-(1/5)*log(500))

thetaz<-(-1441440)*(int$value)
c2k<-1/5
(4*0.06/(thetaz*c2k))^(1/5)*(500)^(-1/5)