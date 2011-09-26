### hierarchical: binomial survival of females, Poisson offspring
# 25 Sep 2011

# parms
x = 3 # when x=8, Pr of survival=0.000; x=2, Pr=0.5; x=-8, Pr=0.999
r = 1

#p = 1/(1 + exp(x))
#R = exp(r)

# simulate
Nt = sort(rep(c(0,1,2,4,8,16),10)) # females introduced
Ot = rbinom(length(Nt), Nt, 1/(1 + exp(x)) ) # females surviving to oviposit
Nt1 = rpois(length(Nt), lambda = exp(r) * Ot) # galls



# specify neg log likelihood function
#parms = list(x = 1.3, r = 1)

nL = function(x, r, debug=FALSE) {
	p = 1/(1 + exp(x))
	R = exp(r)
	L = numeric(length(Nt1))
	for(i in 1:length(Nt1)){
		O = 0:max(Nt[i])
		L[i] = sum(
			dbinom(O, size = Nt[i], p = p) *
			dpois(Nt1[i], lambda= O * R)
		)

	}
	negL = -log(prod(L))
	if(debug)
		cat(x, r, negL, '\n')
	return(negL)
}





# fit
library(bbmle)
#nL(x = 0, r=1)

m1 = mle2(nL, start=list(x=0, r=1), control=list(trace=FALSE))
cat('True values:', 'p =', 1/(1 + exp(x)), 'R =', exp(r), '\n')
cat('Estimates:', 'p =', 1/(1 + exp(coef(m1)['x'])), 'R =', exp(coef(m1)['r']), '\n')
