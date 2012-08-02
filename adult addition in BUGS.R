## Adult addition experiment in BUGS
## 5 April 2012

library(rjags)
library(R2jags)
library(R2WinBUGS)

# loop for power sim
save = data.frame(R.lower=0, R.upper=0, p.lower=0, p.upper=0)

for(i in 1:100){

# exptl set up
F = rep(c(1,2,4,8,16), 12)
nPlants = length(F)

# parms
p = 0.75
R = 3

# simulation
O = rbinom(nPlants, F, p)
G = rpois(nPlants, R * O)

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
		lambda[i] <- R * O[i]
		O[i] ~ dbin(p, F[i])
	}
}
", fill = TRUE)
sink()

# Initial values
inits = function() list(p = 0.5, R = 3.1, O = rep(1, nPlants))

# parms
parms = c('p', 'R')

# data
data = list(F=F, nPlants=nPlants, G=G)

# MCMC controls
ni = 50000
nb = 0.1 * ni


# use rjags
mrj = jags.model('adult-addition.txt', data, inits, n.chains=1,
	n.adapt=1000)

# burn in
update(mrj, 1000)

# take MCMC samples
orj = coda.samples(mrj, parms, n.iter=ni)

save[i,] = matrix(t(HPDinterval(orj)[[1]]), nrow=1)
}

plot(save[,1], ylim=c(min(save[,1:2]), max(save[,1:2])))
points(save[,2], pch=20)
abline(h=R)

plot(save[,3], ylim=c(min(save[,3:4]), max(save[,3:4])))
points(save[,4], pch=20)
abline(h=p)

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











# WinBUGS

# ob = bugs(data, inits, parms, 'adult-addition.txt', 
	# n.thin=1, n.chains=2, n.burnin=nb, n.iter=ni,
	# debug=TRUE)


# JAGS with R2jags

#or2j = jags(data, inits, parms, 'adult-addition.txt', n.chains=2,
	#n.thin=1, n.burnin=nb, n.iter=ni)

#str(or2j$BUGSoutput$sims.list)

#plot(or2j$BUGSoutput$sims.list$R ~ or2j$BUGSoutput$sims.list$p)
#cor(or2j$BUGSoutput$sims.list$R, or2j$BUGSoutput$sims.list$p)
