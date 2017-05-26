rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')
library('Matrix')

setwd('U:\\independent studies\\std geospatial\\github_std_geospatial')
source('gibbs std geospatial functions.R')
source('gibbs std geospatial wrapper.R')
dat=read.csv('fake data.csv',as.is=T)

tipos=c('',2)
nomes.cov.village=c(paste('nl',tipos,sep=''),
        paste('elevation',tipos,sep=''),
        paste('mean_evi',tipos,sep=''),
        paste('dist_water',tipos,sep=''),
        paste('dist_urb',tipos,sep=''))
nomes.cov.ind='agekr.yr'
ngibbs=100

res=gibbs.wrapper(dat,nomes.cov.village,nomes.cov.ind,ngibbs,include.pred=T)
