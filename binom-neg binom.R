### hierarchical binomial-negative binomial
# 26 Sep 2011

# binomal female survival followed by neg binom distribution of offspring

nll.bnb = function(x, r, s, debug=FALSE) {
	p = 1/(1 + exp(x))
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
	if(debug)
		cat(x, r, s, negL, '\n')
	return(negL)
}


BNB = mle2(nL2, start=list(x = 3, r = 1, s=1))
cat('Estimates:', 'p =', 1/(1 + exp(coef(BNB)['x'])), 'R =', exp(coef(BNB)['r']), '\n')


