# ~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/
#	plotting VAA 2010 results.R

###        plotting VAA 2010 results         ###
### density / host-plant quality experiment  ###

# graphical exploration
# plotting number of galls collected from each experimental plant 
# in the 2009-2010 density/quality experiment at Valentine.

# load data, entered gall numbers on map data csv
d = read.csv('~/Documents/DATA/2011 DATA/field/VAA map.csv')

library(ggplot2)
library(rgl)
library(scatterplot3d)
library(bbmle)

# plotting questions:
#		1. should I include the 2 dying plants?
#		2. should I include plants to which I introduced 0 females?


## first 2-dimensional plots
# load someone's hack to allow suppression of top and right borders
#source("~/Documents/R/ggplot2/rnc_ggplot2_border_themes.r")

# 1. plot number of galls in spring 2011 ~ number of females introduced in summer 2010
p1 = ggplot(data = d, aes(x = females, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + 
	scale_x_continuous('females') +
	scale_y_continuous('galls produced') +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank()) # remove gridlines
				
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p1.smooth = p1 + stat_smooth(method = 'lm', colour='royalblue') 
p1.smooth
# save plot as pdf
ggsave('~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/figs/galls~females_lm.pdf', 
	plot = p1.smooth, width=4, height=4)

# 2. plot number of galls in spring 2011 ~ 
#		number of naturally occuring, experimentally removed 2010 galls
p2 = ggplot(data = d, aes(x = natural.galls, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + theme_bw() + 
	scale_x_continuous('previous galls') +
	scale_y_continuous('galls produced') +
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank())
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p2.smooth = p2 + stat_smooth(method = 'lm', colour='royalblue')
p2.smooth
# save plot as a pdf
ggsave(
	'~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/figs/galls~naturalgalls_lm.pdf', 
	plot = p2.smooth, width=4, height=4)

# 3. plot number of galls in spring 2011 ~ time corrected water potential
#		d$psi.r is residuals of linear regression psi ~ time
p3 = ggplot(data = d, aes(x = psi.r, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, size = 2.5) + 
	scale_y_continuous(limits = c(0,12)) +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank())
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p3.smooth = p3 + stat_smooth(method = 'loess', colour='royalblue', span=1)

# save plot as pdf
ggsave('~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/figs/galls~rpsi.pdf', 
	plot = p3.smooth, width=4, height=4)

## second 3-d plots

# 1. plot galls on plants in spring 2011 ~ removed natural galls and introduced females
# highlight.3d colors points in relation to y axis, the number of females
# type='h' displays the height of each point
s3d1 = scatterplot3d(d$natural.galls, d$females, d$galls2011, pch=16, highlight.3d=TRUE, 
	type='h', zlab = 'Galls in year t+1', ylab = 'Introduced females', 
	xlab = 'Galls in year t')
# add a regression plane
fit.pois = mle2(galls2011 ~ dpois(lambda = a + b * natural.galls + c * females), 
	start = list(a = mean(d$galls2011), b = 0, c = 0), method = 'SANN', data=d)
s3d1$plane3d(fit.pois)
dev.print(pdf, file=
	"~/Documents/Analysis repos/Eutreta-diana-density-quality-experiment/figs/3d VAA 2010.pdf",
	width=5.1, height=5.1, pointsize=8) # width/height are in inches


# 2. interactive version of above plot
plot3d(d$galls2011 ~ d$natural.galls + d$females, type='s', size=1, zlab='Number of galls', xlab=
	'Previous year gall density',ylab='Number of females')
plot3d(d$galls2011 ~ d$natural.galls + d$females, type='h', add=TRUE)






######## 2-d: only plants with >0 galls
# 1. plot number of galls in spring 2011 ~ number of females introduced in summer 2010
dplus = d[d$galls2011>0,]

p1 = ggplot(data = dplus, aes(x = females, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + 
	scale_x_continuous('females') +
	scale_y_continuous('galls produced') +
	theme_bw() + 
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank()) # remove gridlines
				
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p1.smooth = p1 + stat_smooth(method = 'lm', colour='royalblue') 
p1.smooth

# 2. plot number of galls in spring 2011 ~ 
#		number of naturally occuring, experimentally removed 2010 galls
p2 = ggplot(data = dplus, aes(x = natural.galls, y = galls2011)) +
	geom_point(colour = 'royalblue', alpha = 1/2, position = position_jitter(w = 0.15, h = 0.15),
	size = 2.5) + theme_bw() + 
	scale_x_continuous('previous galls') +
	scale_y_continuous('galls produced') +
	opts( panel.grid.minor = theme_blank(), panel.grid.major = theme_blank())
		
# now add a smoother with confidence interval, methods include loess, lm, glm (family='poisson')
p2.smooth = p2 + stat_smooth(method = 'lm', colour='royalblue')
p2.smooth





