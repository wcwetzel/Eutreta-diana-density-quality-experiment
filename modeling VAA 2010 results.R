

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

