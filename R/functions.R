



#mid level functions

glasso_opt=function(rho,S.list,full=FALSE,quick=FALSE,
                    start="cold",cov.last=NULL,prec.last=NULL){
  P=dim(S.list[[1]])[1]
  J=length(S.list)
  eps=1e-10
  maxit=10
  array.S=array(unlist(S.list),c(P,P,J))
  S.init=cov2cor(apply(array.S,c(1,2),mean))
  weights=rep(1,J)/J

  if(start=="warm"){
    A=glasso::glasso(S.init,rho=rho,penalize.diagonal=FALSE,thr=1e-8,
                     start = "warm",w.init=cov.last,wi.init=prec.last)
    cov.last=A$w
    prec.last=A$wi
  }else{
    A=glasso::glasso(S.init,rho=rho,penalize.diagonal=FALSE,thr=1e-8)
  }

  
  Sigma.gl=Theta.gl=list()
  Sigma.gl[[1]]=A$w ;Theta.gl[[1]]=A$wi

  if(!(quick)){
    count=1
    diff=eps+1
    while((diff>eps & count<maxit)& any(!is.na(Theta.gl[[count]])))
    {
      weights=plyr::laply(S.list,L.icov.prop,Theta=Theta.gl[[count]])
      weights=weights/sum(weights)
      #print(1/sum(weights^2))
      count=count+1
      Sigma.gl[[count]]=cov2cor(apply(array.S,c(1,2),weighted.mean,w=weights))

      if(start=="warm"){
        A=glasso::glasso(Sigma.gl[[count]],rho=rho,penalize.diagonal=FALSE,thr=1e-8,
                         start = "warm",w.init=cov.last,wi.init=prec.last)
        cov.last=A$w
        prec.last=A$wi
      }else{
        A=glasso::glasso(Sigma.gl[[count]],rho=rho,penalize.diagonal=FALSE,thr=1e-8)
      }
      
      
      Sigma.gl[[count]]=A$w; Theta.gl[[count]]=A$wi
      if(any(!is.na(Theta.gl[[count]]))){
        diff=sum(((Theta.gl[[count]]-Theta.gl[[count-1]])^2)/(P^2))
      } else {diff=sum(((Sigma.gl[[count]]-Sigma.gl[[count-1]])^2)/(P^2))
      }
      #print(diff)
    }

  }
  if(full){
    return(A)
  }else{
    return(A$wi)
  }
}












approx.graph<- function(manyglm.obj,lambdas,calc.lik=FALSE,quick=quick)
{
  res=residuals(manyglm.obj)#this will not work for ordinal data at all
  S=cov2cor(cov0(res))
  P=dim(S)[1]
  N=dim(manyglm.obj$y)[1]
  
  Th.out=list()
  n.lam=length(lambdas)
  A=glasso_opt(lambdas[1],list(S),full=TRUE,quick=FALSE)
  cov.last=A$w
  prec.last=A$wi
  Th.out[[1]]=prec.last
  
  for(i.lam in 2:n.lam){
    A=glasso_opt(lambdas[i.lam],list(S),full=TRUE,quick=FALSE,
               start="warm",cov.last=cov.last,prec.last=prec.last)
    cov.last=A$w
    prec.last=A$wi
    Th.out[[i.lam]]=prec.last
  }
  
  
  
  k=rep(NA,length(Th.out))
  for(j in 1:length(Th.out)){
    k[j]=(sum(Th.out[[j]]!=0)-P)/2
  }
  k.frac=k/(P*(P-1)/2)
  BIC.out=logL=NULL
  if(calc.lik){
    logL=plyr::laply(Th.out,quick.cop.llik,res=res)
    BIC.out=k*log(N)-2*logL-sum(manyglm.obj$two.loglike)
  }
  return(list(k.frac=k.frac,BIC=BIC.out,logL=logL))
}



full.graph.many<- function(manyglm.obj,lambdas,res)
{
  S.list=res$S.list
  res=res$res
  P=dim(S.list[[1]])[1]
  N=dim(manyglm.obj$fitted.values)[1]
  
  Th.out=Sig.out=list()
  n.lam=length(lambdas)
  A=glasso_opt(lambdas[1],S.list,full=TRUE,quick=FALSE)
  cov.last=A$w
  prec.last=A$wi
  Th.out[[1]]=prec.last
  Sig.out[[1]]=cov.last
  
    for(i.lam in 2:n.lam){
      A=glasso_opt(lambdas[i.lam],S.list,full=TRUE,quick=FALSE,
                   start="warm",cov.last=cov.last,prec.last=prec.last)
      cov.last=A$w
      prec.last=A$wi
      Th.out[[i.lam]]=prec.last
      Sig.out[[i.lam]]=cov.last
    }

  
  
  
  k=rep(NA,length(Th.out))
  for(j in 1:length(Th.out)){
    k[j]=(sum(Th.out[[j]]!=0)-P)/2
  }
  k.frac=k/(P*(P-1)/2)
  BIC.out=logL=NULL
  

  logL=plyr::laply(Th.out,ll.icov.all,S.list=S.list,n=N)
  BIC.out=k*log(N)-2*logL#-sum(manyglm.obj$two.loglike)

  
  return(list(k.frac=k.frac,BIC=BIC.out,logL=logL,Th.out=Th.out,Sig.out=Sig.out))
  
}

#low level functions

simulate.res.S<-function(manyglm.obj,n.res=200){
  many.output=list(manyglm.obj)[rep(1,n.res)]
  res=plyr::llply(many.output,residuals)#this will not work for ordinal data at all
  if(min(res[[1]])>-1e-5){#{residuals function for manyany currently output on uniform scale
    res=plyr::llply(res,qnorm)
  }
  S.list=plyr::llply(res,function(x) cov2cor(cov0(x)))#cov0 assumes 0 mean
  return(list(res=res,S.list=S.list))
}

cov0=function(Y){
  n=dim(Y)[1]
  (t(Y)%*%Y)/(n-1)
}








quick.cop.llik<-function(Theta,res){
  P=dim(res)[2]
  N=dim(res)[1]
  
  llik=0
  for(i in 1:N){
    data=matrix(res[i,],P,1)
    llik=llik+1/2*log(det(Theta)) -1/2 *( t(data)%*%( Theta-diag(P) )%*%data )
  }
  
  llik
}


ll.icov.all<-function(theta,S.list,n){
  P=dim(S.list[[1]])[1]
  K=length(S.list)
  inside=0
  for(i in 1:K){
    scatter=(n-1)*S.list[[i]]
    inside=inside+exp(-sum(diag(scatter%*%(theta-diag(P))))/2)
  }
  
  ll=-n/2*log(2*pi) +n/2*log(det(theta)) +log(inside/n)
  #notice theta is sigma inverse hence the first minus is a plus
}


print.saint=function(obj) 
{
  print(paste0("Species interactions from copula model", obj$formula))
  print("Precision matrix")
  print(obj$best.graph$prec)
}

L.icov.prop=function(S,Theta)
{
  out=exp(-1/2*(sum(S*(Theta-diag(dim(Theta)[1])))))#the last bit is the trace 
  return(out)
}


