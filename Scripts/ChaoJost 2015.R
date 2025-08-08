# R scripts for computing diversity (Hill numbers) profile using individual-based abundance data or sampling-unit-based incidence data.
# In all functions, param x is a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (for incidence data).
# For incidence data, the first entry of x must be the number of sampling units.
# In all functions, param q is the diversity order; the suggested range for q is [0, 3].
# If you use the scripts for publishing papers, please cite Chao and Jost 2015 MEE paper (Appendix S8).

#-----------------------------------------------
# Diversity profile estimator (abundance data)
#-----------------------------------------------
#' Chao_Hill_abu(x, q) is a function of obtaining estimators of Hill numbers of order q based on abundance data.
#' @param x a vector of species sample frequencies.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @return a numerical vector of diversity.

Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#-----------------------------------------------
# Diversity profile estimator (incidence data)
#-----------------------------------------------
#' Chao_Hill_inc(x, q) is a function of obtaining estimators of Hill numbers of order q based on incidence data.
#' @param x a vector of species incidence-based sample frequencies. The first entry of x must be the number of sampling units.
#' @param q a numeric or a vector of diversity order.
#' @return a numerical vector.

Chao_Hill_inc = function(x,q){
  n = x[1]
  x = x[-1];x = x[x>0]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/U*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/U*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)*U/n
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,((n/U)^q*A)^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      ((n/U)^q*(A+B))^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#' Chao_Hill(x, q,datatype) combines Chao_Hill_abu and Chao_Hill_inc given a specified datatype (either abundance data or incidence data).
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Chao_Hill = function(x,q,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    est = Chao_Hill_abu(x,q)
  }else{
    est = Chao_Hill_inc(x,q)
  }
  return(est)
}

#-----------------------
# The empirical profile
#-----------------------
#' Hill(x, q, datatype) is a function of obtaining the empirical Hill numbers of order q based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Hill <- function(x,q,datatype = c("abundance","incidence")){
  if(datatype=="incidence"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

#-----------------------------------------
# The bootstrap method for obtaining s.e.
#-----------------------------------------
#' Bt_prob_abu(x) is a function of estimating the species probabilities in the bootstrap assemblage based on abundance data.
#' @param x a vector of species sample frequencies.
#' @return a numeric vector.

Bt_prob_abu = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  W = (1-C)/sum(x/n*(1-x/n)^n)

  p.new = x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0 = (1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob_inc(x) is a function of estimating the species incidence probabilities in the bootstrap assemblage based on incidence data.
#' @param x a vector of incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @return a numeric vector.

Bt_prob_inc = function(x){
  n = x[1]
  x = x[-1]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  C=1-f1/U*(1-A)
  W=U/n*(1-C)/sum(x/n*(1-x/n)^n)

  p.new=x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0=U/n*(1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob(x,datatype) combines the two functions Bt_prob_abu and Bt_prob_inc for a specified datatype.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numeric vector.

Bt_prob = function(x,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    prob = Bt_prob_abu(x)
  }else{
    prob = Bt_prob_inc(x)
  }
  return(prob)
}

#' Bootstrap.CI(x,q,B,datatype,conf) is a function of calculating the bootsrapping standard error based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param B an integer to specify the number of replications in the bootstrap procedure, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200)..
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 3 matrices including respectively the difference between the average and lower confidence bound of the B bootstrap estimates,
#'  the difference between the upper confidence bound and the average of the B bootstrap estimates, and the bootstrap standard error of the diversity estimate.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

Bootstrap.CI = function(x,q,B = 1000,datatype = c("abundance","incidence"),conf = 0.95){
  datatype = match.arg(datatype,c("abundance","incidence"))
  p.new = Bt_prob(x,datatype)
  n = ifelse(datatype=="abundance",sum(x),x[1])
  # set.seed(456)
  if(datatype=="abundance"){
    data.bt = rmultinom(B,n,p.new)
  }else{
    data.bt = rbinom(length(p.new)*B,n,p.new)
    data.bt = matrix(data.bt,ncol=B)
    data.bt = rbind(rep(n,B),data.bt)
  }

  mle = apply(data.bt,2,function(x)Hill(x,q,datatype))
  pro = apply(data.bt,2,function(x)Chao_Hill(x,q,datatype))

  mle.mean = rowMeans(mle)
  pro.mean = rowMeans(pro)

  LCI.mle =  -apply(mle,1,function(x)quantile(x,probs = (1-conf)/2)) + mle.mean
  UCI.mle = apply(mle,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - mle.mean

  LCI.pro =  -apply(pro,1,function(x)quantile(x,probs = (1-conf)/2)) + pro.mean
  UCI.pro = apply(pro,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - pro.mean

  LCI = rbind(LCI.mle,LCI.pro)
  UCI = rbind(UCI.mle,UCI.pro)

  sd.mle = apply(mle,1,sd)
  sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
  se = rbind(sd.mle,sd.pro)

  return(list(LCI=LCI,UCI=UCI,se=se))

}

#------------------------
# Main function ChaoHill
#------------------------
#' ChaoHill(dat, datatype, from, to, interval, B, conf) is the function of calculating the empirical and the proposed diversity profile,
#' their bootsrap standard errors and confidance intervals.
#' @param dat a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param from a numeric number of diversity order q (the start order of profile).
#' @param to a numeric number of diversity order q (the end order of profile).
#' @param interval a numeric number to specify each increment of q from the start to end order.
#' @param B an integer to specify the number of bootstrap replications, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200).
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 4 matrices including respectively diversity estimates, bootstrap standard errors, lower confidence bounds, and upper confidence bounds.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

ChaoHill <- function(dat, datatype=c("abundance", "incidence"), from=0, to=3, interval=0.1, B=1000, conf=0.95){
  datatype = match.arg(datatype,c("abundance","incidence"))
  # for real data estimation

  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  q <- seq(from, to, by=interval)

  #-------------
  #Estimation
  #-------------
  MLE=Hill(dat,q,datatype)

  qD_pro=Chao_Hill(dat,q,datatype)

  CI_bound = Bootstrap.CI(dat,q,B,datatype,conf)
  se = CI_bound$se
  #-------------------
  #Confidence interval
  #-------------------
  tab.est=data.frame(rbind(MLE,qD_pro))

  LCI <- tab.est - CI_bound$LCI
  UCI <- tab.est + CI_bound$UCI

  colnames(tab.est) <- colnames(se) <- colnames(LCI) <- colnames(UCI) <- paste("q = ", q, sep="")
  rownames(tab.est) <- rownames(se) <- rownames(LCI) <- rownames(UCI) <- c("Observed", "Chao_2013")
  return(list(EST = tab.est,
              SD = se,
              LCI = LCI,
              UCI = UCI))

}


#----------------------------
# Plot of confidence interval
#----------------------------
#' conf.reg(x_axis,LCL,UCL,...) is a function to plot the confidence region.
#'
#' \code{conf.reg} uses polygon to draw a confidence band plot
#'
#' @param x_axis a vector of diversity orders.
#' @param LCL a vector of lower confidence bounds.
#' @param UCL a vector of upper confidence bounds.
#' @param ... further arguments to be passed to \code{polygon}
#' @return a polygon plot

conf.reg=function(x_axis,LCL,UCL,...) {
  x.sort <- order(x_axis)
  x <- x_axis[x.sort]
  LCL <- LCL[x.sort]
  UCL <- UCL[x.sort]
  polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
}
