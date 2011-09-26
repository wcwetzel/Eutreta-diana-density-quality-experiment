# ~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/
#	modeling VAA 2010 results.R

###        modeling VAA 2010 results         ###
### density / host-plant quality experiment  ###

library(bbmle)
library(emdbook)

# load data, entered gall numbers on map data csv
d = read.csv('~/Documents/DATA/2011 DATA/field/VAA map.csv')

d$area = d$d1 * d$d2
d$volume = d$area * d$h
d$gpf = d$galls2011 / d$females
d$gpf[is.nan(d$gpf)] = 0
d$lgpf = log(d$gpf)

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

# m05. females, no density dependence
m05 = mle2(galls2011 ~ dnbinom( mu = r * females, size=s ), 
	start = list(r = 0.5, s=0.4), data = d)

# m075 females, galls, no DD
m075 = mle2(galls2011 ~ dnbinom( mu =  females * exp(r * natural.galls), size=s ), 
	start = list(r = 0.5, s=0.4, b = 0), data = d)

# m1. females, ricker -- it can't seem to calculate profiles for this parameterization
#m1 = mle2(galls2011 ~ dnbinom( mu = females * exp(r * (1 - females / k)), size=s ), 
#	start = list(r = -0.3, k = 4, s=0.4), data = d)
	
# m1. females, ricker
m1 = mle2(galls2011 ~ dnbinom( mu = R * females * exp(-b * females), size=s ), 
	start = list(R = 0.5, b = 1, s=0.4), data = d)
	
# m2. females, beverton holt
m2a = mle2(galls2011 ~ dnbinom( mu = exp(r) * females / (1 + a * females), size=s ), 
	start = list(r = 0.8, a = 0.1, s=0.4), data = d)
	
# m2. females, beverton holt
m2b = mle2(galls2011 ~ dnbinom( mu = exp(r) * females / (1 + exp(a) * females), size=s ), 
	start = list(r = 0.8, a = 0.1, s=0.4), data = d)
	
# m2. females, beverton holt
m2c = mle2(galls2011 ~ dnbinom( mu = exp(r) * females / (1 + females), size=s ), 
	start = list(r = 0.8, s=0.4), data = d)
	
# m2. females, beverton holt
m2d = mle2(galls2011 ~ dnbinom( mu = exp(r) * females / (a + females), size=s ), 
	start = list(r = 0.8, a = 1, s=0.4), data = d)

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




AICtab(m0, m05, m075, m1, m2a, m2b, m2c, m2d, m3, m4, m5)
BICtab(m0, m05, m075, m1, m2a, m2b, m2c, m2d, m3, m4, m5)




## compounded binomial - negative binomial ##

















#### PLOTTING PREDICTIONS AND CIs ####
# adding fits and CIs to plots
# first run 'plotting VAA 2010 results.R'


# mean predictions and CI for m0.5
# predicted values
fp = 0:18
pred.m05 = data.frame(females = fp, galls = fp * coef(m05)['r'])

# confidence intervals
# nbinom is easy: mu is a parm
# for m05, mu = r * females 
# so the only uncertainty in mu comes from r
profile.m05 = profile(m05)
ci.m05 = confint(profile.m05)
pred.m05$ymin = ci.m05['r', 1] * fp
pred.m05$ymax = ci.m05['r', 2] * fp

p1.fitted = p1 + geom_smooth(aes(x=fp, y=galls, ymin=ymin, ymax=ymax), 
	data=pred.m05, stat='identity', colour='royalblue', fill='royalblue')

# mean predictions and CI for m1
pred.m1 = data.frame(females = fp, galls = coef(m1)['R'] * fp * exp(-coef(m1)['b'] * fp ))

# confidence intervals for m1
# for m1, mu = females * exp(r * (1 - females / k))
# which values of r and k do I use for upper and lower bounds?
# I imagine, the bounds that give the lowest low and highest high?
profile.m1 = profile(m1)
ci.m1 = confint(profile.m1)
pred.m1$ymin = ci.m1['R', 1] * fp * exp(-ci.m1['b', 2] * fp)
pred.m1$ymax = ci.m1['R', 2] * fp * exp(0 * fp)

p1.fitted = p1.fitted + geom_smooth(aes(x=fp, y=galls, ymin=ymin, ymax=ymax), 
	data=pred.m1, stat='identity', colour='red', fill='red', alpha=1/4)

print(p1.fitted)

# mean predictions and CI for m2
pred.m2c = data.frame(females = fp, galls = coef(m2c)['r'] * fp / (1 + fp))

# confidence intervals for m2
# for m2, mu = r * females / (1 + a * females)
# which values of r and a do I use for upper and lower bounds?
# I imagine, the bounds that give the lowest low and highest high?
profile.m2c = profile(m2c)
ci.m2c = confint(profile.m2c)
pred.m2c$ymin = ci.m2c['r', 1] * fp / (1 + fp) 
pred.m2c$ymax = ci.m2c['r', 2] * fp / (1 + fp)

p1.fitted = p1 + geom_smooth(aes(x=fp, y=galls, ymin=ymin, ymax=ymax), 
	data=pred.m2c, stat='identity', colour='yellow', fill='yellow', alpha=1/4)

print(p1.fitted)



gp.m2 = fp * coef(m2)['r'] / (1 + coef(m2)['a'] * fp)

p1.smooth2 = p1.smooth + geom_smooth(aes(x=fp, y = gp.m1), colour='red') +
	geom_smooth(aes(x=fp, y=gp.m2), colour='green3', lty=2) +

print(p1.smooth2)










################### OLD ###################################
mf1 = mle2(d$galls2011 ~ dpois(lambda = a * d$females), data=d, start=list(a = 1/4))
plot(galls2011 ~ females, data=d)
abline(mf1)
mf2 = mle2(d$galls2011 ~ dnbinom(mu = a * d$females, size = s), 
	start=list(a = 1/4, s = 1), data=d)
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

