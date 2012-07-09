# ~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/
#	modeling VAA 2010 results.R

###        modeling VAA 2010 results         ###
### density / host-plant quality experiment  ###

library(bbmle)
#library(emdbook)
library(rethinking)

# load data, entered gall numbers on map data csv
d = read.csv('~/Documents/DATA/2011 DATA/field/VAA map.csv')

d$area = d$d1 * d$d2
d$volume = d$area * d$h /1000
d$gpf = d$galls2011 / d$females
d$gpf[is.nan(d$gpf)] = 0
d$lgpf = log(d$gpf)
d$I[d$females==0] = 0
d$I[d$females>0] = 1


#------- do galls increase with females or natural galls? --------#
# 7 May 2012 # 
# only with cages with >0 females
df = d[d$females>0,]


# females #
plot(galls2011 ~ females, data=df)

# intercept
m0 = mle2(galls2011 ~ dnbinom(mu = mu, size=s),
	start=list(mu = 2.5, s=1), data=df)

# no intercept linear females (no density dependence)
m1 = mle2(galls2011 ~ dnbinom( mu = R * females, size=s ), 
	start = list(R = 0.5, s=0.4), data = df)

# intercept linear females (no density dependence)
m2 = mle2(galls2011 ~ dnbinom( mu = mu + R * females, size=s ), 
	start = list(mu = 2.5, R = 0.5, s=0.4), data = df)

abline(h=coef(m0)['mu'])
abline(a=0, b=coef(m1)['R'])
abline(a=coef(m2)['mu'], b=coef(m2)['R'])


AICctab(m0, m1, m2, nobs=nrow(df))
anova(m0,m2)


# natural galls #
# intercept
m0 = mle2(galls2011 ~ dnbinom(mu = mu, size=s),
	start=list(mu = 2.5, s=1), data=df)

# intercept linear females (no density dependence)
m1 = mle2(galls2011 ~ dnbinom( mu = mu + R * natural.galls, size=s ), 
	start = list(mu = 2.5, R = 0.5, s=0.4), data = df)

plot(galls2011 ~ natural.galls, data=df)
abline(h=coef(m0)['mu'])
abline(a=coef(m1)['mu'], b=coef(m1)['R'])


AICctab(m0, m1, nobs=nrow(df))
anova(m0,m2)


#-------------------end 7 may 2012------------------#





# using approach from qualifying exam
# see 'qe presentation.pdf'
# Ricker between number of females and number of galls (larvae)
# quality variables added in exponent
# N(t+1) = females * exp(r(1-females/k) + b1*quality1 + b2*quality2)
# negbinom or zero-inflated nbinom stochasticity
# going with nbinom instead of pois b/c var is several times the mean
# also it makes sense that the lambda on each shrub would be gamma distributed


# m00. intercept model, makes much more sense than plain intercept
# 0 galls if 0 females, mu galls if >0 females
m00 = mle2(galls2011 ~ dnbinom(mu = mu * I, size=s),
	start=list(mu = 2.5, s=1), data=d)

# m10. linear females (no density dependence)
m10 = mle2(galls2011 ~ dnbinom( mu = R * females, size=s ), 
	start = list(R = 0.5, s=0.4), data = d)

confint(m10)
postm10 = sample.naive.posterior(m10)
new.females = 0:18
mu = sapply(new.females, function(z) mean(postm10[,1] * z))
mu.ci = sapply(new.females, function(z) HPDI(postm10[,1] * z))

plot(galls2011 ~ females, data=d)
lines(new.females, mu)
lines(new.females, mu.ci[1,], lty=2)
lines(new.females, mu.ci[2,], lty=2)

# m20. nonlinear (ricker) females
m20 = mle2(galls2011 ~ dnbinom( mu = R * females * exp(-k * females), size=s ), 
	start = list(R = 0.5, k = 1, s=0.4), data = d)

postm20 = sample.naive.posterior(m20)
mu20 = sapply(new.females, function(z) mean(postm10[,1] * z * exp(-postm10[,2] * z)))
mu20.ci = sapply(new.females, function(z) HPDI(postm10[,1] * z * exp(-postm10[,2] * z)))

plot(galls2011 ~ females, data=d)
lines(new.females, mu20, col='steelblue')
lines(new.females, mu20.ci[1,], lty=2, col='blue')
lines(new.females, mu20.ci[2,], lty=2, col='blue')

# m01. intercept females, log linear natural galls
m01 = mle2(galls2011 ~ dnbinom(mu = exp(a + b * natural.galls) * I, size=s),
	start=list(a = 2.5, b=0, s=1), data=d)

# m11. linear females, log linear natural galls
m11 = mle2(galls2011 ~ dnbinom( mu = exp(a + b * natural.galls) * females, size=s ), 
	start = list(a = 0.5, b=0, s=0.4), data = d)

# m21. nonlinear (ricker) females, log linear natural galls
m21 = mle2(galls2011 ~ dnbinom( mu = exp(a + b * natural.galls) * 
	females * exp(-k * females), size=s ), start = list(a = 0.5, b = 1, k = 1, s=0.4), data = d)

# m02. intercept females, log linear psi.r
m02 = mle2(galls2011 ~ dnbinom(mu = exp(a + b * psi.r) * I, size=s),
	start=list(a = 2.5, b=0, s=1), data=d)

# m12. linear females, log linear psi.r
m12 = mle2(galls2011 ~ dnbinom( mu = exp(a + b * psi.r) * females, size=s ), 
	start = list(a = 1, b = 0, s=0.4), data = d)

# m22. nonlinear (ricker) females, log linear psi.r
m22 = mle2(galls2011 ~ dnbinom( mu = exp(a + b * psi.r) * females * exp(-k * females), size=s ), 
	start = list(a = 0.5, b = 1, k = 1, s=0.4), data = d)
	
# m1b. linear females, log linear natural galls, log linear psi.r
m1b = mle2(galls2011 ~ dnbinom( mu = exp(a + b * natural.galls + c * psi.r) * females, size=s ), 
	start = list(a = 0.5, b=0, c=0, s=0.4), data = d)

# m00size. intercept females, log linear natural galls
m00size = mle2(galls2011 ~ dnbinom(mu = exp(a + b * volume) * I, size=s),
	start=list(a = 2.5, b=0, s=0.5), data=d)
	
# m10size. linear females (no density dependence)
m10size = mle2(galls2011 ~ dnbinom( mu = exp(a + b * volume) * females, size=s ), 
	start = list(a = 0.5, b=0, s=0.4), data = d)

AICctab(m00, m10, m20, m01, m11, m21, m02, 
	m12, m22, m1b, m00size, m10size, weights=TRUE)
BICtab(m00, m10, m20, m01, m11, m21, m02, 
	m12, m22, m1b, m00size, m10size, weights=TRUE, nobs=30)


## compounded binomial-nbinom ##

# m00 with bnb. intercept females >0
nll.m00bnb = function(mu, x, s) {
	p = 1/(1 + exp(x)) # this transformation is necessary for profiling binom-nbinom
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = mu * d$I, size = s)
		)
	}
	negL = -log(prod(L))	
	cat('mu=', mu, 'x=', x, 'p=', 1/(1+exp(x)), 's=', s, '-nll=', negL, '\n')
	return(negL)
}
m00bnb = mle2(nll.m00bnb, start=list(mu=2, x=0, s=0.5), data=d)
profile.m00bnb = profile(m00bnb)
plot(profile.m00bnb)


# as above but with Poisson
nll.m00bp = function(mu, x) {
	p = 1/(1 + exp(x)) # this transformation is necessary for profiling binom-nbinom
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dpois(d$galls2011[i], lambda = mu * d$I)
		)
	}
	negL = -log(prod(L))	
	cat('mu=', mu, 'x=', x, 'p=', 1/(1+exp(x)), '-nll=', negL, '\n')
	return(negL)
}
m00bp = mle2(nll.m00bp, start=list(mu=2, x=0), data=d)
profile.m00bp = profile(m00bp)
confint(profile.m00bp)

# m10 with bnb, linear females
nll.m10bnb = function(R, x, s) {
	p = 1/(1 + exp(x)) # logistic transformation for survival probability
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = R * females, size = s)
		)
	}
	negL = -log(prod(L))
}
m10bnb = mle2(nll.m10bnb, start=list(R=1, x=0, s=0.5), data=d)
profile.m10bnb = profile(m10bnb)

# m20 with bnb, nonlinear (ricker) females
nll.m20bnb = function(r, k, x, s, debug=TRUE) {
	p = 1/(1 + exp(x)) # logistic transformation for survival probability
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = exp(r) * females * exp(-k * females), size = exp(s))
		)
	}
	negL = -log(prod(L))
	if(debug) 
		cat('r:', r, 'R:', log(r), 'k:', k, 'x:', x, 'p:', 1/(1+exp(x)), 's:', s, negL, '\n')
	return(negL)
}
m20bnb = mle2(nll.m20bnb, start=list(r=0, k=1, x=0, s=0.5), data=d)

# m01 with bnb. intercept females >0, linear natural galls
nll.m01bnb = function(a, b, x, s) {
	p = 1/(1 + exp(x)) # this transformation is necessary for profiling binom-nbinom
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = exp(a + b * natural.galls) * d$I, size = s)
		)
	}
	negL = -log(prod(L))
}
m01bnb = mle2(nll.m01bnb, start=list(a=2, b=0, x=0, s=0.5), data=d)
profile.m01bnb = profile(m01bnb)


nll.m11bnb = function(a, b, x, s) {
	p = 1/(1 + exp(x)) # logistic transformation for survival probability
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = exp(a + b * natural.galls) * females, size = s)
		)
	}
	negL = -log(prod(L))
}
m11bnb = mle2(nll.m11bnb, start=list(a=1, b=0, x=0, s=0.5), data=d)




AICtab(m00, m10, m20, m01, m11, m21, m02, m12, m22, m1b, 
	m00bnb, m10bnb, m20bnb, m01bnb, m11bnb, weights=TRUE, nobs=30)
BICtab(m00, m10, m20, m01, m11, m21, m02, 
	m12, m22, m1b, BNB, weights=TRUE, nobs=30)




# old BNB with debugging code
NLL.BNB = function(x, r, s, debug=FALSE) {
	p = 1/(1 + exp(x)) # these transformations are necessary for profiling binom-nbinom
	R = exp(r)
	L = numeric(length(d$galls2011))
	for(i in 1:length(L)){
		O = 0:max(d$females[i])
		L[i] = sum(
			dbinom(O, size = d$females[i], p = p) *
			dnbinom(d$galls2011[i], mu = O * R, size = s)
		)
	}
	negL = -log(prod(L))
	if(debug) cat(x, r, s, negL, '\n')
	return(negL)
}


BNB = mle2(NLL.BNB, start=list(x = 1, r = 1, s=1))
#cat('Estimates:', 'p =', 1/(1 + exp(coef(BNB)['x'])), 'R =', exp(coef(BNB)['r']), '\n')

profile.BNB = profile(BNB)
ci.BNB = confint(profile.BNB)
1/ (1 + exp(ci.BNB[1,]))
exp(ci.BNB[2,])
plot(profile.BNB)

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

ggsave(filename = 
	'~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/figs/galls~females_fits.pdf',
	plot = p1.fitted)

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





