# ~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/
#	modeling VAA 2010 results - more than zero.R

###        modeling VAA 2010 results         ###
### density / host-plant quality experiment  ###
###       only plants with > 0 galls         ###

library(bbmle)
#library(emdbook)

# load data, entered gall numbers on map data csv
d = read.csv('~/Documents/DATA/2011 DATA/field/VAA map.csv')
d = d[d$galls2011>0,]

d$area = d$d1 * d$d2
d$volume = d$area * d$h /1000
d$gpf = d$galls2011 / d$females
d$gpf[is.nan(d$gpf)] = 0
d$lgpf = log(d$gpf)
d$I[d$females==0] = 0
d$I[d$females>0] = 1

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
m00 = mle2(galls2011 ~ dnbinom(mu = mu, size=s),
	start=list(mu = 2.5, s=1), data=d)

# m10. linear females (no density dependence)
m10 = mle2(galls2011 ~ dnbinom( mu = R * females, size=s ), 
	start = list(R = 0.5, s=0.4), data = d)

# m20. nonlinear (ricker) females
m20 = mle2(galls2011 ~ dnbinom( mu = R * females * exp(-k * females), size=s ), 
	start = list(R = 0.5, k = 1, s=0.4), data = d)

# m01. intercept females, log linear natural galls
m01 = mle2(galls2011 ~ dnbinom(mu = exp(a + b * natural.galls), size=s),
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
	m12, m22, m1b, m00size, m10size)
BICtab(m00, m10, m20, m01, m11, m21, m02, 
	m12, m22, m1b, m00size, m10size, weights=TRUE, nobs=30)