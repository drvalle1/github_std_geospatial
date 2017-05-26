compare1=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango)
}

seq1=1000:4000

compare1(alpha.true,colMeans(vec.alpha[seq1,]))
compare1(gamma.true,colMeans(vec.gamma[seq1,]))

#look at gamma intercept
plot(vec.gamma[seq1,1],type='l')
abline(h=gamma.true[1],col='red')

plot(density(vec.gamma[seq1,1]),type='l')
abline(v=gamma.true[1],col='red')

#look at rho: 
plot(vec.outros[seq1,'rho'],type='l')
abline(h=rho.true)

plot(density(vec.outros[seq1,'rho'],from=0))
abline(v=rho.true)

acf(vec.outros[seq1,'rho'])

#look at sig2
plot(vec.outros[seq1,'sig2'],type='l')
abline(h=sigma2.true,col='red')

plot(density(vec.outros[seq1,'sig2'],from=0))
abline(v=sigma2.true)

acf(vec.outros[seq1,'sig2'])

cor(vec.outros,use="pairwise.complete.obs")

#look at beta
plot(vec.betas[seq1],type='l')
abline(h=betas.true,col='red')
