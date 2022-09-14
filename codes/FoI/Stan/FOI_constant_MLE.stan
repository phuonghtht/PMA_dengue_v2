data {

	int<lower=0> a;// number of age groups, split by year
 
  int datap[6,a]; //vector of sero status (negative, primary (DENV1-4), and secondary for all age groups), split by year

	int veca[a]; //vector of age constitute of specific age according to year rather than ages in order.
}

parameters {

real <lower=0, upper=1> lambda;

}

transformed parameters {
//the probabilities of being neg, primary and secondary at each age group 

real pneg[a];// the probabilities of neg at each age group and year
real pprimsero1[a]; // the probabilities of primary infection of each serotype at each age group and year
real pprimsero2[a];
real pprimsero3[a];
real pprimsero4[a];
real psec[a];// the probabilities of each secondary infection at each age group and year

// store log-likelihood and vector of probabilities for each category
real log_like[a] ;
vector[6] vecProb ;

// calculating all these probabilities
for (ar in 1:a) {

pneg[ar]=exp(-lambda*veca[ar]*4);// the prob of neg at each age group and year

pprimsero1[ar]=(exp(-3*lambda*veca[ar]))*(1-exp(-lambda*veca[ar]));
pprimsero2[ar]=pprimsero1[ar];
pprimsero3[ar]=pprimsero1[ar];
pprimsero4[ar]=pprimsero1[ar];

psec[ar]= 1 - pneg[ar] - pprimsero1[ar] - pprimsero2[ar] - pprimsero3[ar] - pprimsero4[ar];
// the probabilities of each secondary infection at each age group and year

vecProb[1] = pneg[ar] ;
vecProb[2] = pprimsero1[ar] ;
vecProb[3] = pprimsero2[ar] ;
vecProb[4] = pprimsero3[ar] ;
vecProb[5] = pprimsero4[ar] ;
vecProb[6] = psec[ar] ;

log_like[ar] = multinomial_lpmf(datap[1:6,ar] | vecProb[1:6]);
}
}

model {
//vector to put the model results in
vector[6] vp;


//prior on the parameter we are estimating
lambda~ beta(2,5); 

// the estimates of the things we have data for given the proportions above.

for (ar in 1:a){
vp[1]=pneg[ar];
vp[2]=pprimsero1[ar];
vp[3]=pprimsero2[ar];
vp[4]=pprimsero3[ar];
vp[5]=pprimsero4[ar];
vp[6]=psec[ar];
//likelihood thing assume multinomial
datap[1:6, ar]  ~  multinomial(vp) ;
}

}

generated quantities {
	real sumloglike;
	sumloglike = 0;	
	for (ar in 1:a) {
		sumloglike += log_like[ar];
		}
}

