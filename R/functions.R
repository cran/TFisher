#' @importFrom stats pnorm qnorm uniroot

########################
#Skew Normal Estimation#
########################
findsnalpha = function(skewness){
  skewness = sign(skewness)*min(0.9952717,abs(skewness))
  delta = sqrt(pi/2*(skewness^2)^(1/3)/((skewness^2)^(1/3)+(2-pi/2)^(2/3)))
  return(sign(skewness)*delta/sqrt(1-delta^2))
}

findsnomega = function(variance, alpha){
  delta2 = alpha^2/(1+alpha^2)
  return(sqrt(variance/(1-2*delta2/pi)))
}

findsnxi = function(m, alpha, omega){
  delta = alpha/sqrt(1+alpha^2)
  return(m-omega*delta*sqrt(2/pi))
}

##################
#Gaussian mixture#
##################
F0 = function(q){
  return(pnorm(q,0,1))
}

F1 = function(q,eps,mu){
  return(eps*pnorm(q,mu,1)+(1-eps)*pnorm(q,0,1))
}

F0_1 = function(p){
  return(qnorm(p,0,1))
}

F1_1 = function(p,eps,mu,br=c(-1000,1000)){
  G = function(q) F1(q,eps,mu) - p
  return(uniroot(G,br)$root)
}

vecF1_1 = Vectorize(F1_1, vectorize.args = "p")

D = function(q, eps, mu){
  return(1-F1(F0_1(1-q), eps, mu))
}

D_1 = function(q, eps, mu){
  return(1-F0(vecF1_1(1-q, eps, mu)))
}
