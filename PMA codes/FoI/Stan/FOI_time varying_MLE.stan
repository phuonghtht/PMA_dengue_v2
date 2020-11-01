data {
  
  // HCMC infections
  int<lower=0> a;// number of age groups
  int step_year;
  int<lower=0> l_vecj;
  
  int veca[a]; //vector age constitute of specific age rather than ages in order.
  
  int vecy[a]; //vector year of collecting samples
  
  int vecj[l_vecj];// vector vecj comprises integer value => for index of lambda
  
  int datap[3,a]; ////for 2013 this will be three groups (primary, secondary and negative for all age groups)
  
  
  
}

parameters {
  
  vector<lower=0,upper=1>[l_vecj/step_year] lambda; // 2017-2013+max(age of 2013): 34/2=17 years HC, KH: 35/2=18 yrs
  
}

transformed parameters {
  //the probabilities of being neg, primary and secondary at each age group 
  ////// HCMC
  // need to specify what type of data and how these are...
  
  real pneg[a];// the probabilities of neg at each age group and year
  real pprimsero[a] ; // the probabilities of primary infection of each serotype at each age group and year
  real psec[a];// the probabilities of each secondary infection at each age group and year
  real exposure;
  
  // store log-likelihood and vector of probabilities for each category
  real log_like[a] ;
  vector[3] vecProb ;  
  // calculating all these probabilities
  for ( ar in 1:a){
    int startLambda=2017-vecy[ar]+1; // IdxLambda start from 2017 is 1,2,3,....
    int birthyear=vecy[ar]-veca[ar];
    int endLambda=2017-birthyear;

    exposure = sum(lambda[vecj[startLambda:endLambda]]); // expectation of any dengue exposure
    pneg[ar] = exp(-exposure*4);// the prob of neg at each age group and year
    pprimsero[ar] = 4*(1-exp(-exposure))*exp(-exposure*3);
    psec[ar]= 1-pneg[ar]-pprimsero[ar];
    
    vecProb[1] = pneg[ar] ;
    vecProb[2] = pprimsero[ar] ;
    vecProb[3] = psec[ar] ;
    
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




