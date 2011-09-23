###           VAA 2010 results               ###
### density / host-plant quality experiment  ###

# 1st plotting and then model fitting for number of galls collected from each
# experimental plant in the 2009-2010 density/quality experiment at Valentine.

# load data, entered gall numbers on map dataframe
d = read.csv('/Users/will/Documents/DATA/2011 DATA/field/VAA map.csv')

library(ggplot2)
library(lattice)
library(rgl)
library(fields)
library(scatterplot3d)
library(bbmle)


## first 2-dimensional plots
# load someone's hack to allow suppression of top and right borders
source("/Users/will/Documents/R/ggplot2/rnc_ggplot2_border_themes.r")

# 1. plot number of galls in spring 2011 ~ number of females introduced in summer 2010
p1 = ggplot(data = d, aes(x = females, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
		size = 2.5) +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(),
		panel.border = theme_border(c('left', 'bottom')))
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p1.smooth = p1 + stat_smooth(method = 'loess', colour='royalblue', span=1) 
print(p1.smooth) # display the plot

# 2. plot number of galls in spring 2011 ~ 
#		number of naturally occuring, experimentally removed 2010 galls
p2 = ggplot(data = d, aes(x = natural.galls, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
		size = 2.5) +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(),
		panel.border = theme_border(c('left', 'bottom')))
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p2.smooth = p2 + stat_smooth(method = 'loess', colour='royalblue', span=1)
print(p2.smooth) # display the plot


plot(galls2011[d$females>0 & d$shrubdying2011==0] ~ natural.galls[d$females>0 & d$shrubdying2011==0],
	data=d)
plot(galls2011[d$shrubdying2011==0] ~ females[d$shrubdying2011==0], data=d)

#cloud( galls2011 ~ natural.galls + females, data=d)

plot3d(d$galls2011 ~ d$natural.galls + d$females, type='s', size=1, zlab='Number of galls', xlab=
	'Previous year gall density',ylab='Number of females')
plot3d(d$galls2011 ~ d$natural.galls + d$females, type='h', add=TRUE)

plot3d(d$galls2011[d$females>0 & d$shrubdying2011==0] ~ d$natural.galls[d$females>0 & d$shrubdying2011==0]
	+ d$females[d$females>0 & d$shrubdying2011==0], type='s', size=1, )
plot3d(d$galls2011[d$females>0 & d$shrubdying2011==0] ~ d$natural.galls[d$females>0 & d$shrubdying2011==0] 
	+ d$females[d$females>0 & d$shrubdying2011==0], type='h', add=TRUE)

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

