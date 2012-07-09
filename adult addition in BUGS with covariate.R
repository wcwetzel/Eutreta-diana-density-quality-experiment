# Adult addition experient in BUGS with covariate
# 5 April 2012


## Adult addition experiment in BUGS
## 5 April 2012

library(rjags)

save = data.frame(R.lower=0, R.upper=0, p.lower=0, p.upper=0)

for(i in 1:100){
	
# exptl set up
F = rep(c(1,2,4,8,16), 12)
Q = sort(rep(1:5, 12))
nPlants = length(F)

# parms
p = 0.5
R = 2

# simulation
O = rbinom(nPlants, F, p)
G = rpois(nPlants, R * Q * O)

# Specify model in BUGS language
sink('adult-addition.txt')
cat("
model {
	# Priors
	p ~ dunif(0, 1)
	R ~ dunif(0,10)
	
	# Likelihood
	for(i in 1:nPlants){
		G[i] ~ dpois(lambda[i])
		lambda[i] <- R * O[i] * Q[i]
		O[i] ~ dbin(p, F[i])
	}
}
", fill = TRUE)
sink()

# Initial values
inits = function() list(p = 0.5, R = 1.6, O = rep(1, nPlants))

# parms
parms = c('p', 'R')

# data
data = list(F=F, nPlants=nPlants, G=G, Q=Q)

# MCMC controls
ni = 50000
nb = 0.1 * ni

# use rjags #

# compile model
mrjc = jags.model('adult-addition.txt', data, inits, n.chains=1,
	n.adapt=1000)

# burn in
update(mrjc, 1000)

# take MCMC samples
orj = coda.samples(mrjc, parms, n.iter=ni)

save[i,] = matrix(t(HPDinterval(orj)[[1]]), nrow=1)
}

plot(save[,1], ylim=c(min(save[,1:2]), max(save[,1:2])))
points(save[,2], pch=20)
abline(h=R)

plot(save[,3], ylim=c(min(save[,3:4]), max(save[,3:4])))
points(save[,4], pch=20)
abline(h=p)



# analyze samples
plot(orj, ask=TRUE)
summary(orj)
HPDinterval(orj)
gelman.diag(orj)
autocorr.plot(orj)
crosscorr(orj)
crosscorr.plot(orj)
acfplot(orj)
densplot(orj)
effectiveSize(orj)
gelman.plot(orj)
rejectionRate(orj)

orjmat = as.data.frame(as.matrix(orj))
plot(orjmat$R ~ orjmat$p)