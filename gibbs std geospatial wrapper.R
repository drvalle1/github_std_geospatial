gibbs.wrapper=function(dat,nomes.cov.village,nomes.cov.ind,ngibbs,include.pred){
  dat=dat[order(dat$loc.id),]
  
  #for prediction
  ind.pred=which(is.na(dat$microsc1))
  npred=length(ind.pred)
  
  nomes=c('LATNUM','LONGNUM','loc.id',nomes.cov.village)
  clust=unique(round(dat[,sort(nomes)],4))
  clust=clust[order(clust$loc.id),]
  nclust=nrow(clust); nclust
  
  #get dmat
  dmat=create.dmat(dat)
  dmat=Matrix(dmat)
  dtd=t(dmat)%*%dmat
  
  #get wmat and gamma
  ind=which(colnames(clust)%in%c('LATNUM','LONGNUM','loc.id'))
  wmat=data.matrix(cbind(1,clust[,-ind]))
  gamma=rep(0,ncol(wmat))
  
  #get xmat and betas
  xmat=data.matrix(dat[,nomes.cov.ind])
  betas=rep(0,ncol(xmat))
  xtx=t(xmat)%*%xmat
  
  z=rep(NA,nrow(dat))
  cond=!is.na(dat$microsc1) & dat$microsc1==1; z[cond]=1
  cond=!is.na(dat$microsc1) & dat$microsc1==0; z[cond]=-1
  cond=is.na(dat$microsc1); z[cond]=-0.5
  
  #distance matrix
  dist1=as.matrix(dist(clust[,c('LATNUM','LONGNUM')]))
  nvil=nrow(dist1)
  
  #discretize rho to improve mixing
  cor1=seq(from=0.05,to=0.95,by=0.01)
  rho.vals=-quantile(dist1,0.1)/log(cor1)
  # exp(-quantile(dist1,0.1)/rho.vals)
  # round(exp(-quantile(dist1,0.9)/rho.vals),2)
  ldet.vals=rep(NA,length(cor1))
  Sigma.vals=invSigma.vals=list()
  for (i in 1:length(rho.vals)){
    Sigma.vals[[i]]=exp(-dist1/rho.vals[i])
    invSigma.vals[[i]]=solve(Sigma.vals[[i]])
    ldet.vals[i]=determinant(Sigma.vals[[i]],logarithm = T)$modulus[[1]]
  }
  
  #get initial values
  ind.rho=3
  rho=rho.vals[ind.rho]
  alpha=rmvnorm(1,mean=rep(0,ncol(dist1)),sigma=Sigma.vals[[ind.rho]])
  sig2=2
  
  #useful things for gibbs sampler
  vec.alpha=matrix(NA,ngibbs,length(alpha))
  vec.outros=matrix(NA,ngibbs,2) #rho+sig2
  vec.gamma=matrix(NA,ngibbs,length(gamma))
  vec.betas=matrix(NA,ngibbs,length(betas))
  vec.pred=matrix(NA,ngibbs,npred)
  
  param=list(alpha=t(alpha),rho=rho,gamma=gamma,betas=betas,
             invSigma=invSigma.vals[[ind.rho]],sig2=sig2,z=z)
  
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    param$alpha=update.alpha(param,dtd,xmat,dmat,wmat) #t(alpha.true) #
    
    tmp=update.rho(param,wmat,rho.vals,invSigma.vals,nvil,ldet.vals)
    param$rho=tmp$rho
    param$invSigma=tmp$invSigma
    
    param$sig2=update.sig2(param,nvil,wmat) #sigma2.true 
    param$gamma=update.gamma(param,wmat)#gamma.true
    
    param$betas=update.betas(param,xtx,xmat,dmat)#betas.true
    param$z=update.z(param,xmat,dat) #z.true#
    
    #store results
    vec.alpha[i,]=param$alpha
    vec.outros[i,]=c(param$rho,param$sig2) 
    vec.gamma[i,]=param$gamma
    vec.betas[i,]=param$betas
    
    if (include.pred){
      media=param$alpha[dat$loc.id]+xmat%*%param$betas
      vec.pred[i,]=pnorm(media[ind.pred])
    }
  }
  
  colnames(vec.outros)=c('rho','sig2')
  colnames(vec.gamma)=colnames(wmat)
  colnames(vec.betas)=colnames(xmat)  
  
  list(alpha=vec.alpha,outros=vec.outros,gamma=vec.gamma,betas=vec.betas,pred=vec.pred)
}
