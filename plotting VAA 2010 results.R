# /Users/will/Documents/Analysis repos/
# Eutreta-diana-density-quality-experiment/plotting VAA 2010 results.R

###        plotting VAA 2010 results         ###
### density / host-plant quality experiment  ###

# graphical exploration
# plotting number of galls collected from each experimental plant 
# in the 2009-2010 density/quality experiment at Valentine.

# load data, entered gall numbers on map dataframe
d = read.csv('/Users/will/Documents/DATA/2011 DATA/field/VAA map.csv')

library(ggplot2)
library(lattice)
library(rgl)
library(fields)
library(scatterplot3d)
library(bbmle)
library(emdbook)

# plotting questions:
#		1. should I include the 2 dying plants?
#		2. should I include plants to which I introduced 0 females?


## first 2-dimensional plots
# load someone's hack to allow suppression of top and right borders
source("/Users/will/Documents/R/ggplot2/rnc_ggplot2_border_themes.r")

# 1. plot number of galls in spring 2011 ~ number of females introduced in summer 2010
p1 = ggplot(data = d, aes(x = females, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + 
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(), # remove gridlines
		panel.border = theme_border(c('left', 'bottom'))) # remove top and right of plot box
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p1.smooth = p1 + stat_smooth(method = 'loess', colour='royalblue', span=1) 
pdf('/Users/will/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/galls~females.pdf')
print(p1.smooth)
dev.off()

# 2. plot number of galls in spring 2011 ~ 
#		number of naturally occuring, experimentally removed 2010 galls
p2 = ggplot(data = d, aes(x = natural.galls, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(),
		panel.border = theme_border(c('left', 'bottom')))
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p2.smooth = p2 + stat_smooth(method = 'loess', colour='royalblue', span=1)

# save plot as a pdf
ggsave('/Users/will/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/galls~naturalgalls2.pdf', plot = p2.smooth, width=4, height=4)

# 3. plot number of galls in spring 2011 ~ time corrected water potential
#		d$psi.r is residuals of linear regression psi ~ time
p3 = ggplot(data = d, aes(x = psi.r, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, size = 2.5) + 
	scale_y_continuous(limits = c(0,12)) +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(),
		panel.border = theme_border(c('left', 'bottom')))
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p3.smooth = p3 + stat_smooth(method = 'loess', colour='royalblue', span=1)
pdf('/Users/will/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/galls~rpsi.pdf')
print(p3.smooth)
dev.off()

## second 3-d plots

# 1. plot galls on plants in spring 2011 ~ removed natural galls and introduced females
# highlight.3d colors points in relation to y axis, the number of females
# type='h' displays the height of each point
pdf('/Users/will/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/3d-VAA2010.pdf')
s3d1 = scatterplot3d(d$natural.galls, d$females, d$galls2011, pch=16, highlight.3d=TRUE, 
	type='h', zlab = 'Galls in year t+1', ylab = 'Introduced females', 
	xlab = 'Galls in year t')
# add a regression plane
fit.pois = mle2(galls2011 ~ dpois(lambda = a + b * natural.galls + c * females), 
	start = list(a = mean(d$galls2011), b = 0, c = 0), method = 'SANN', data=d)
s3d1$plane3d(fit.pois)
dev.off()

# 2. interactive version of above plot
plot3d(d$galls2011 ~ d$natural.galls + d$females, type='s', size=1, zlab='Number of galls', xlab=
	'Previous year gall density',ylab='Number of females')
plot3d(d$galls2011 ~ d$natural.galls + d$females, type='h', add=TRUE)











