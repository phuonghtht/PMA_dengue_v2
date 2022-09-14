data {
  int<lower=0> a;// number of age groups by years
  int step_year;
  int<lower=0> l_vecj;
  int veca[a]; //vector age constitute of specific age, useful in case missing age.
  int vecy[a]; //vector year, specific year of data is indicated.
  //int lamb_y; // number of years that FoI is estimated in the model
  //int n_y;  // number of years that FoI can be estimated from the data.
  int vecj[l_vecj]; // vector vecj comprises integer value => for index of lambda in case of grouping
    
  int datap[6,a]; ////for 2013 this will be sixgroups (neg,primary(1,2,3,4),and secondary  for all age groups)
}

parameters {
// CHECKING LAMBDA VALUES
//  real <lower=0, upper=1> lambda;
  vector <lower=0,upper=1>[l_vecj/step_year] lambda; // 2017-1984 =34/5=7years group
}

transformed parameters {
  real pneg[a];// the probabilities of neg at each age group and year
  real pprimsero1[a]; // the probabilities of primary infection of each serotype at each age group and year
  real pprimsero2[a];
  real pprimsero3[a];
  real pprimsero4[a];
  real psec[a];// the probabilities of each secondary infection at each age group and year
  real exposure;

  // store log-likelihood and vector of probabilities for each category
  real log_like[a] ;
  vector[6] vecProb ;   
    
  // calculating all these probabilities
    
  for (ar in 1:a) {
    int startLambda=2017-vecy[ar]+1; // IdxLambda start from 2017 is 1,2,3,....
    int birthyear=vecy[ar]-veca[ar];
    int endLambda=2017-birthyear+1;
    
// CHECKING LAMBDA VALUES
//pneg[ar]=exp(-lambda*veca[ar]*4);// the prob of neg at each age group and year
//pprimsero1[ar]=(exp(-3*lambda*veca[ar]))*(1-exp(-lambda*veca[ar]));
//pprimsero2[ar]=pprimsero1[ar];
//pprimsero3[ar]=pprimsero1[ar];
//pprimsero4[ar]=pprimsero1[ar];    
    
    exposure = sum(lambda[vecj[startLambda:endLambda]]); // expectation of any dengue exposure
    pneg[ar] = exp(-exposure*4);// the prob of neg at each age group and year
    //pprimsero1[ar] = 4*(1-exp(-exposure))*exp(-exposure*3);
    pprimsero1[ar] = (1-exp(-exposure))*exp(-exposure*3);
    pprimsero2[ar] = pprimsero1[ar];
    pprimsero3[ar] = pprimsero1[ar];
    pprimsero4[ar] = pprimsero1[ar];
    psec[ar]= 1 - pneg[ar] - pprimsero1[ar] - pprimsero2[ar] - pprimsero3[ar] - pprimsero4[ar];
        
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

// CHECKING LAMBDA VALUES
//lambda ~ beta(2,5);

  for (i in 1:(l_vecj/step_year)) {
    lambda[i] ~ beta(2,5); // lambda is a matrix => vector distribution
  }

  // the estimates of the things we have data for given the proportions above.

  for (ar in 1:a) {
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

