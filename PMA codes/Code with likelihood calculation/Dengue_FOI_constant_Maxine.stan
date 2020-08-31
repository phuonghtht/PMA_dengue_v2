data {

	int<lower=0> a;// number of age groups
 
	int datap[3,a]; //vector of sero status (negative, primary, and secondary for all age groups)

	int veca[a]; //vector of age constitute of specific age rather than ages in order.

	int vectotal[a]; //vector of total data points per age group in dataset
}

parameters {

real<lower=0, upper=1> lambda;

}

transformed parameters {
//the probabilities of being neg, primary and secondary at each age group 

real pneg[a];// the probabilities of neg at each age group and year
real pprimsero[a] ; // the probabilities of primary infection of each serotype at each age group and year
real psec[a];// the probabilities of each secondary infection at each age group and year
real log_like[a];
vector[3] vecProb;

// calculating all these probabilities
for (ar in 1:a) {

pneg[ar]=exp(-lambda*veca[ar]*4);// the prob of neg at each age group and year

pprimsero[ar]=4*exp(-lambda*veca[ar]*3)*(1-exp(-lambda*veca[ar])); // the probabilities of primary infection of each serotype at each age group and year//primmary infection regardless of numbers of time of infection with such serotype.

psec[ar]= 1-pneg[ar]-pprimsero[ar];// the probabilities of each secondary infection at each age group and year

vecProb[1] = pneg[ar];
vecProb[2] = pprimsero[ar];
vecProb[3] = psec[ar];

log_like[ar] = multinomial_lpmf(datap[1:3,ar] | vecProb[1:3]);
}
}

model {
//vector to put the model results in
vector[3] vp;


//prior on the parameter we are estimating
lambda~ beta(2,5); 

// the estimates of the things we have data for given the proportions above.

for (ar in 1:a){
vp[1]=pneg[ar];
vp[2]=pprimsero[ar];
vp[3]=psec[ar];
//likelihood thing assume multinomial
datap[1:3, ar]  ~  multinomial(vp) ;
}

}

generated quantities {
	real sumloglike;
	sumloglike = 0;	
	for (ar in 1:a) {
		sumloglike += log_like[ar];
		}
}

