rm(list=ls(all=TRUE))
set.seed(123)
library('mvtnorm')

#get data
setwd('U:\\ghana\\fever and malaria\\single DHS survey')
dat=read.csv('actual data.csv',as.is=T)

#create village level matrix
tipos=c('',2)
nomes=c('LATNUM','LONGNUM','loc.id',
        paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
clust=unique(round(dat[,sort(nomes)],4))
clust=clust[order(clust$loc.id),]; nrow(clust)

#create distance matrix
dist1=as.matrix(dist(clust[,c('LATNUM','LONGNUM')]))
hist(dist1)
rho.true=rho=2
seq1=seq(from=0,to=max(dist1),length.out=1000)
plot(seq1,exp(-seq1/rho),type='l')
abline(v=1:6,col='grey')
exp(-median(dist1)/rho)

#generate alphas
ind=which(colnames(clust)%in%c('LATNUM','LONGNUM','loc.id'))
env1=clust[,-ind]
Sigma=exp(-dist1/rho)
sigma2.true=sigma2=2
wmat=cbind(1,data.matrix(env1))
gamma=c(-0.3,rnorm(ncol(wmat)-1,mean=0,sd=0.5))
gamma.true=gamma
media=wmat%*%gamma
range(media)
alpha.true=alpha=rmvnorm(1,media,sigma2*Sigma)

#generate infection status
dat1=dat
xmat=data.matrix(dat1[,'agekr.yr'])
betas.true=betas=runif(ncol(xmat),min=-0.4,max=0.4)
media=alpha[dat1$loc.id]+xmat%*%betas
range(xmat%*%betas)
z.true=z=rnorm(nrow(dat1),mean=media,sd=1)

#replace infected
mean(dat1$microsc1,na.rm=T)
dat1$microsc1=ifelse(z.true>0,1,0)
mean(dat1$microsc1,na.rm=T)

setwd('U:\\independent studies\\std geospatial\\github_std_geospatial')
write.csv(dat1,'fake data.csv',row.names=F)
