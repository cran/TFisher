#' @importFrom stats pnorm qnorm dnorm uniroot pgamma
#' @importFrom mvtnorm pmvnorm

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

###############################
#Variance-Covariance Estimation
###############################
getTFisherMean = function(n, tau1, tau2){
  2*n*tau1*(1-log(tau1/tau2))
}

getTFisherCovM = function(M, TAU1, TAU2){
  K = length(TAU1)
  #outer(1:K,1:K,function(r,c)getTFisherCov_v1(M, TAU1[r], TAU2[r], TAU1[c], TAU2[c]))
  COV = matrix(NA, K, K)
  for(i in 1:K){
    for(j in i:K){
      COV[i,j] = getTFisherCov_v1(M, TAU1[i], TAU2[i], TAU1[j], TAU2[j])
    }
  }
  COV[lower.tri(COV)] = t(COV)[lower.tri(COV)]
  return(COV)
}

getTFisherCov_v1 = function(M, tau11, tau21, tau12, tau22){
  #calculate cov(Ui, Uj), where Ui=Xi^2I(Xi^2>b), Uj=Xj^2I(Xj^2>b)
  #where cor(Xi,Xj)=rho_ij, b=qnorm(1-tau1/2)
  #variance = sum_{ij}cov(Ui, Uj)
  n = length(M[1,])

  #Define Constants for W1
  b_1 = qnorm(1-tau11/2)
  b2_1 = b_1^2
  phib_1 = dnorm(b_1)
  bphib_1 = b_1*phib_1
  m0_1 = 1-tau11
  m2_1 = m0_1 - 2*bphib_1
  m4_1 = 3*m2_1 - 2*b2_1*bphib_1
  varU_1 = 3-m4_1 - (1-m2_1)^2 #VAR=var(Ui)=var(Uj)
  varY_1 = 4*tau11*(1+(1-tau11)*(1-log(tau11/tau21))^2)

  #Define Constants for W2
  b_2 = qnorm(1-tau12/2)
  b2_2 = b_2^2
  phib_2 = dnorm(b_2)
  bphib_2 = b_2*phib_2
  m0_2 = 1-tau12
  m2_2 = m0_2 - 2*bphib_2
  m4_2 = 3*m2_2 - 2*b2_2*bphib_2
  varU_2 = 3-m4_2 - (1-m2_2)^2 #VAR=var(Ui)=var(Uj)
  varY_2 = 4*tau12*(1+(1-tau12)*(1-log(tau12/tau22))^2)

  varRatio = sqrt(varY_1/varU_1*varY_2/varU_2)

  # Integration Domain
  x = seq(-8,8,length.out=1000)
  dx = x[2] - x[1]
  x2 = x^2; x1trunc = x2*(x2>b2_1) #x2trunc = x2*(x2>b2_2)#
  phix = dnorm(x)

  # Extract duplicate rho
  tbl = table(M[col(M)!=row(M)]) # off-diagonals
  RHO = as.numeric(names(tbl))
  FREQ = as.numeric(tbl)
  tau1min = min(tau11,tau12)

  COV = 4*n*tau1min*(1+(1+log(tau21/tau1min))*(1+log(tau22/tau1min))-tau11*tau12/tau1min*(1+log(tau21/tau11))*(1+log(tau22/tau12)))
  if(n==1){
    return(COV)
  }else{
    for(i in 1:length(RHO)){
      #################################
      rho = RHO[i]                    #
      rho2 = rho^2                    #
      sqrho2 = sqrt(1-rho2)           #
      xx = -rho*x                     #
      A = xx/sqrho2; B = b_2/sqrho2     #
      fx = A+B; fx2 = fx^2; fx3 = fx^3#
      gx = A-B; gx3 = gx^3            #
      #first part above 20% time#######
      #################################
      phigx = dnorm(gx)               #
      phifx = dnorm(fx)               #
      Phigx = pnorm(gx)               #
      Phifx = pnorm(fx)               #
      #second part above 75% time######
      #################################
      m0x = Phifx - Phigx             #
      m1x = phigx - phifx             #
      m2x = m0x + m1x*fx - 2*B*phigx  #
      m2x = m0x + gx*phigx - fx*phifx #
      #third part above 5% time########

      hx = phix*(rho2*x2*m0x + (1-rho2)*m2x + 2*rho*sqrho2*x*m1x)
      covU = 2*rho2*b2_1*bphib_1 + (1+2*rho2-(tau12+2*bphib_2))*(tau11+2*bphib_1) - sum(x1trunc*hx*dx)
      covY = covU*varRatio
      COV = COV + FREQ[i]*covY
    }
    return(COV)
  }
}

getTFisherVar_v1 = function(M, tau1, tau2){
  #calculate cov(Ui, Uj), where Ui=Xi^2I(Xi^2>b), Uj=Xj^2I(Xj^2>b)
  #where cor(Xi,Xj)=rho_ij, b=qnorm(1-tau1/2)
  #variance = sum_{ij}cov(Ui, Uj)
  n = length(M[1,])

  #Define Constants
  b = qnorm(1-tau1/2)
  b2 = b^2
  phib = dnorm(b)
  bphib = b*phib
  m0 = 1-tau1
  m2 = m0 - 2*bphib
  m4 = 3*m2 - 2*b2*bphib
  varU = 3-m4 - (1-m2)^2 #VAR=var(Ui)=var(Uj)
  varY = 4*tau1*(1+(1-tau1)*(1-log(tau1/tau2))^2)
  varRatio = varY/varU

  # Integration Domain
  x = seq(-8,8,length.out=1000)
  dx = x[2] - x[1]
  x2 = x^2; x2trunc = x2*(x2>b2)
  phix = dnorm(x)

  # Extract duplicate rho
  tbl = table(M[col(M)!=row(M)]) # off-diagonals
  RHO = as.numeric(names(tbl))
  FREQ = as.numeric(tbl)

  if(n==1){
    return(varY)
  }else{
    VAR = n*varY # the variance under independence
    for(i in 1:length(RHO)){
      #################################
      rho = RHO[i]                    #
      rho2 = rho^2                    #
      sqrho2 = sqrt(1-rho2)           #
      xx = -rho*x                     #
      A = xx/sqrho2; B = b/sqrho2     #
      fx = A+B; fx2 = fx^2; fx3 = fx^3#
      gx = A-B; gx3 = gx^3            #
      #first part above 20% time#######
      #################################
      phigx = dnorm(gx)               #
      phifx = dnorm(fx)               #
      Phigx = pnorm(gx)               #
      Phifx = pnorm(fx)               #
      #second part above 75% time######
      #################################
      m0x = Phifx - Phigx             #
      m1x = phigx - phifx             #
      m2x = m0x + m1x*fx - 2*B*phigx  #
      m2x = m0x + gx*phigx - fx*phifx #
      #third part above 5% time########

      hx = phix*(rho2*x2*m0x + (1-rho2)*m2x + 2*rho*sqrho2*x*m1x)
      covU = 2*rho2*b2*bphib + (1+2*rho2-(tau1+2*bphib))*(tau1+2*bphib) - sum(x2trunc*hx*dx)
      covY = covU*varRatio
      VAR = VAR + FREQ[i]*covY
    }
    return(VAR)
  }
}


q.tfisher_N <- function(p, n, tau1, tau2, M){
  n = length(M[1,])
  #truncation parameter
  bound = qnorm(1-tau1/2)
  p0 = pmvnorm(lower=-rep(bound,n),upper=rep(bound,n),sigma=M)[1]
  if(p<p0){
    return(0)
  }else{
    mu = getTFisherMean(n, tau1, tau2)
    sigma2 = getTFisherVar_v1(M, tau1, tau2)

    Fw = function(q)pnorm((q-mu)/sqrt(sigma2))
    lower = 0
    upper = max(2*n*(1+log(tau2/tau1))+6*sqrt(4*n*tau1*(1+(1-tau1)*(1+log(tau2/tau1))^2)))
    maxnumIter = 1000
    numIter = 1
    p_cal_lower = Fw(lower)
    p_cal_upper = Fw(upper)

    while((p_cal_lower - p) * (p_cal_upper - p) > 0 && numIter <
          maxnumIter) {
      if(p_cal_lower > p) {
        lower = lower * 1.1
        p_cal_lower = Fw(lower)
      }
      if(p_cal_upper < p) {
        upper = upper * 1.1
        p_cal_upper = Fw(upper)
      }
      numIter = numIter + 1
    }
    if(p_cal_lower>p){
      return(lower)
    }else{
      q = mean(c(lower, upper))
      p_cal = Fw(q)
      error = (p_cal - p)/p
      numIter = 1
      while (abs(error) > 1e-08 && numIter < maxnumIter) {
        if(error > 0) {
          upper = q
          q = mean(c(lower, upper))
          p_cal = Fw(q)
          error = (p_cal - p)/p
        }else{
          lower = q
          q = mean(c(lower, upper))
          p_cal = Fw(q)
          error = (p_cal - p)/p
        }
        numIter = numIter + 1
      }
      if(numIter < maxnumIter) {
        return(q)
      }else{
        return(paste0("The quantile is larger than ", upper))
      }
    }
  }
}
