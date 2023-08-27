#### ARL evaluation via simulation for the new EWMA-type chart for Poisson counts
#######################
ds<-1.0;
lam0<-8.0;
lam<-ds*lam0
gamx<-1
gamy<-57
LCL<-0
UCL<-8
listRL<-vector()
gamxy<-gamx+gamy
it<-5000
for(i in 1:it){
  r<-0
  y<-floor(lam0)
  j<-0
  listy<-c(y)
  while(y>=LCL&UCL>=y){
    x<-rpois(1,lambda=(lam));
    a<-(gamx)*x+gamy*y+r;
    y<-floor(a/gamxy);
    r<-a-gamxy*y;
    j<-j+1;
    }
  
  listRL[i]<-j
  
}
ARLsim<-mean(listRL)
SDRLsim<-sd(listRL)
se<-SDRLsim/sqrt(it)
cat(" ARL:",ARLsim," SDRL:",SDRLsim," se:",se,"\n")
############ END ###################
