### This code evaluates the performance of the new upper-sided EWMA-type chart for Poisson observations
################################################
dewmapARLOut<-function(gamx,gamy,K,lam0,ds){
# gamx and gamy are integer weights of the scheme
# K is the upper control limit
# lam0 is the IC Poisson mean
# ds is the shift in lam0
#########################
# the OoC lambda
  lam<-ds*lam0
# defining a_max and b_max
  amax<-(gamx+gamy)*(K+1)-1
  bmax<-gamx+(K+1)*gamy-1
# Constructing the transition probabilities matrix Q
  Q<-matrix(0,nrow=(bmax+1),ncol=(bmax+1))
  for(b in 0:bmax){
    b1<-b+1
    xmax<-floor((amax-b)/gamx)
    #Q[b1,b2]<-0
    for(x in 0:xmax){
      a<-gamx*x+b
      y<-floor(a/(gamx+gamy))
      #r<-a-(gamx+gamy)*y
      b2<-a-gamx*y+1
      Q[b1,b2]<-Q[b1,b2]+dpois(x,lambda=lam,log=FALSE)
    }
    
  }
# Identity Matrix  
  ID<-diag(x=1,nrow=(bmax+1),ncol=(bmax+1))
# a vector with 1s
  l1<-rep(1,bmax+1)
# definining the initial probabilities vector
  i0<-gamy*floor(lam0)+1
  q<-rep(0,bmax+1)
  q[i0]<-1
# The inverse of the (ID-Q) matrix
  W<-solve(ID-Q)
# ARL
  arl<-q%*%W%*%l1
  arl
}
###################
# A function of the SDRL of the new upper-sided EWMA-type chart for Poisson observations
dewmapSDRLOut<-function(gamx,gamy,K,lam0,ds){
  lam<-ds*lam0
  amax<-(gamx+gamy)*(K+1)-1
  amax
  bmax<-gamx+(K+1)*gamy-1
  bmax
  Q<-matrix(0,nrow=(bmax+1),ncol=(bmax+1))
  for(b in 0:bmax){
    b1<-b+1
    xmax<-floor((amax-b)/gamx)
    for(x in 0:xmax){
      a<-gamx*x+b
      y<-floor(a/(gamx+gamy))
      #r<-a-(gamx+gamy)*y
      b2<-a-gamx*y+1
      Q[b1,b2]<-Q[b1,b2]+dpois(x,lambda=lam,log=FALSE)
    }
    
  }
  ID<-diag(x=1,nrow=(bmax+1),ncol=(bmax+1))
  l1<-rep(1,bmax+1)
  i0<-gamy*floor(lam0)+1
  q<-rep(0,bmax+1)
  q[i0]<-1
  W<-solve(ID-Q)
  nu1<-q%*%W%*%l1
  nu2<-2*q%*%W%*%W%*%l1
# SDRL
  sdrl<-sqrt(nu2-nu1^2+nu1)
sdrl
}
############################
## An example: Evaluating the ARL and SDRL for upward shifts in lambda
## vector of shifts
ds1<-c(1.0,1.1,1.2,1.3,1.5,2.0)
## The OoC ARL and SDRL as function of the shift ds
ARLout<-function(ds){dewmapARLOut(gamx,gamy,K,lam0,ds)}
SDRLout<-function(ds){dewmapSDRLOut(gamx,gamy,K,lam0,ds)}
### 
for(ds in ds1){
ARL<-dewmapOut(gamx=1,gamy=57,K=8,lam0=8,ds)
SDRL<-dewmapSDRLOut(gamx=1,gamy=57,K=8,lam0=8,ds)
  cat(" shift:",ds," ARL:",ARL, " SDRL:",SDRL,"\n")
  #cat(sprintf("shift=%3.1f  arl=%5.2f  sdrl=%5.2f\n",ds,ARLout(ds),SDRLout(ds)))
}
###################### END ############################