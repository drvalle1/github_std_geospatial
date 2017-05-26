rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')
library('Matrix')

# setwd('U:\\ghana\\fever and malaria\\single DHS survey')
# dat=read.csv('actual data.csv',as.is=T)

setwd('U:\\independent studies\\std geospatial\\github_std_geospatial')
source('gibbs std geospatial functions.R')
source('gibbs std geospatial wrapper.R')

tipos=c('',2)
nomes.cov.village=c(paste('nl',tipos,sep=''),
                    paste('elevation',tipos,sep=''),
                    paste('mean_evi',tipos,sep=''),
                    paste('dist_water',tipos,sep=''),
                    paste('dist_urb',tipos,sep=''))
nomes.cov.ind='agekr.yr'
ngibbs=50000
seq1=(9*ngibbs/10):ngibbs

for (valid in 1:10){

  setwd('U:\\independent studies\\gaussian process\\validation')
  nome=paste('valid',valid,'.csv',sep='')
  dat=read.csv(nome,as.is=T)
  ind.pred=which(is.na(dat$microsc1))
  
  res=gibbs.wrapper(dat,nomes.cov.village,nomes.cov.ind,ngibbs,include.pred=T)
  
  setwd('U:\\independent studies\\gaussian process\\validation\\results')
  fim=data.frame(ind=dat$ind[ind.pred],prev.geospat=colMeans(res$pred[seq1,]))
  nome=paste('res',valid,' geospat.csv',sep='')
  write.csv(fim,nome,row.names=F)
}
