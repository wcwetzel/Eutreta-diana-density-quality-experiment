# ~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/
#	modeling VAA 2010 results.R

###        modeling VAA 2010 results         ###
### density / host-plant quality experiment  ###

library(bbmle)
library(emdbook)

# load data, entered gall numbers on map data csv
d = read.csv('~/Documents/DATA/2011 DATA/field/VAA map.csv')

# using approach from qualifying exam
# see 'qe presentation.pdf'
# Ricker between number of females and number of galls (larvae)
# quality variables added in exponent
# N(t+1) = females * exp(r(1-females/k) + b1*quality1 + b2*quality2)
# negbinom or zero-inflated nbinom stochasticity
# going with nbinom instead of pois b/c var is several times the mean
# also it makes sense that the lambda on each shrub would be gamma distributed

# m0. intercept model
m0 = mle2(galls2011 ~ dnbinom( mu = mu, size = s), start = list(mu = 
	mean(d$galls2011), s=1), data=d)

# m0.5. females, no density dependence
m0.5 = mle2(galls2011 ~ dnbinom( mu = r * females, size=s ), 
	start = list(r = 0.5, s=0.4), data = d)

# m0.75 females, galls, no DD
m0.75 = mle2(galls2011 ~ dnbinom( mu = r * females + b * natural.galls, size=s ), 
	start = list(r = 0.5, s=0.4, b = 0), data = d)

# m1. females, ricker
m1 = mle2(galls2011 ~ dnbinom( mu = females * exp(r * (1 - females / k)), size=s ), 
	start = list(r = -0.3, k = 4, s=0.4), data = d)

# m2. females, beverton holt
m2 = mle2(galls2011 ~ dnbinom( mu = r * females / (1 + a * females), size=s ), 
	start = list(r = 0.8, a = 0.1, s=0.4), data = d)

# m3. females, psi.r, ricker
m3 = mle2(galls2011 ~ dnbinom( mu = females * exp(r * (1 - females / k) + 
	b1 * psi.r), size=s ), start = list(r = -0.3, k = 4, s=0.4, b1=0), data = d)

# m4. females, galls, ricker
m4 = mle2(galls2011 ~ dnbinom( mu = females * exp(r * (1 - females / k) + 
	b1 * natural.galls), size=s ), start = list(r = -0.3, k = 4, s=0.4, b1=0), 
	data = d)

# m5. females, galls, psi.r, ricker
m5 = mle2(galls2011 ~ dnbinom( mu = females * exp(r * (1 - females / k) + 
	b1 * natural.galls + b2 * psi.r), size=s ), start = list(r = -0.3, k = 4, s=0.4, 
	b1=0, b2=0), data = d)




AICtab(m0, m0.5, m1, m2, m3, m4, m5)

plot(galls2011 ~ females, data=d)
fp = 0:18
gp = fp * exp(coef(m1)['r'] * (1 - fp / coef(m1)['k']))
points(gp ~ fp, type='l')
p1.smooth2 = p1.smooth + geom_smooth(aes(ymin = , ymax = ))














################### OLD ###################################
mf1 = mle2(d$galls2011 ~ dpois(lambda = a * d$females), data=d, start=list(a = 1/4))
plot(galls2011 ~ females, data=d)
abline(mf1)
mf2 = mle2(d$galls2011 ~ dnbinom(mu = a * d$females, size = s), start=list(a = 1/4, s = 1), data=d)
abline(mf2, col='blue')
mf3 = mle2(d$galls2011 ~ dnbinom(mu = a * d$females / (b + d$females), size=s), 
	start=list(a = 7, b = 2, s=1), data=d, method='SANN')
points(x=1:16, y = coef(mf3)[1] * 1:16 / (coef(mf3)[2] + 1:16), type='l', col='red')

mf4 = mle2(d$galls2011 ~ dnbinom(mu = a * (1 - exp(-b * d$females)), size=s),
	start=list(a = 7, b = 1/4, s = 1), data=d)
points(x=1:16, y = coef(mf4)[1] * (1 - exp(-coef(mf4)[2] * 1:16)), type='l', col='orange')

mf5 = mle2(d$galls2011 ~ dnbinom(mu = exp(b * d$females), size=s), 
	start=list(b = 1, s=1), data=d)
points(x=1:16, y = exp(coef(mf5)[1] * 1:16), type='l', col='green')

mfg1 = mle2(d$galls2011 ~ dpois(lambda = a * (d$females + d$natural.galls)), data=d, 
	start=list(a = 1/4))
abline(mfg1, lty=2)

mfg2 = mle2(d$galls2011 ~ dnbinom(mu = a * (d$females + d$natural.galls), size = s), 
	start=list(a = 1/4, s = 1), data=d)
abline(mf2, col='blue', lty=2)

mg1 = mle2(d$galls2011 ~ dpois(lambda = a * (d$natural.galls+1)), data=d, 
	start=list(a = 1/4))
abline(mfg1, lty=2)

mg2 = mle2(d$galls2011 ~ dnbinom(mu = a * (1 + d$natural.galls), size = s), 
	start=list(a = 1/4, s = 1), data=d)
abline(mf2, col='blue', lty=2)

mf2.1 = mle2(d$galls2011 ~ dnbinom(mu = a * d$females + b * (psi.r+.3), size = s), 
	start=list(a = 1/4, b=1, s = 1), data=d)
abline(mf2, col='blue')

AICtab(mf1, mf2, mf3, mf4, mf5, mfg1, mfg2, mg1, mg2)



# try a zero inflated or other compound model

# zero-inflated nbinom from bolker book package

mz1 = mle2(galls2011 ~ dzinbinom(mu = a * females, size = s, zprob = zprob), 
	start = list(a = 1/4, s = 1, zprob = 0.5), data = d)


AICtab(mf1, mf2, mf3, mf4, mf5, mfg1, mfg2, mg1, mg2, mz1)
##############################################################

