data {
    int<lower=0> a;// number of age groups by years
    int step_year;
    int<lower=0> l_vecj;
    int veca[a]; //vector age constitute of specific age, useful in case missing age.
    int vecy[a]; //vector year, specific year of data is indicated. 
    
    //int lamb_y; // number of years that FoI is estimated in the model ( vary in case of grouping years)
    //int n_y;  // number of years that FoI can be estimated from the data.
    int vecj[l_vecj]; // vector vecj comprises integer value => for index of lambda in case of grouping
    int datap[6,a]; ////for 2013 this will be sixgroups (neg,primary(1,2,3,4),and secondary  for all age groups)
}

parameters {
//  real <lower=0, upper=1> lambda;
  matrix<lower=0,upper=1>[4,l_vecj/step_year] lambda;//
}

transformed parameters {
  real pneg[a];
  real pprimsero1[a] ;
  real pprimsero2[a] ;
  real pprimsero3[a] ;
  real pprimsero4[a] ;
  real psec[a];
    
  // store log-likelihood and vector of probabilities for each category
  real log_like[a] ;
  vector[6] vecProb ;

  real exposure1;
  real exposure2;
  real exposure3;
  real exposure4;
    
  // calculating all these probabilities
    
  for (ar in 1:a){
    int startLambda=2017-vecy[ar]+1; // IdxLambda start from 2017 is 1,2,3,....
    int birthyear=vecy[ar]-veca[ar];
    int endLambda=2017-birthyear+1;

// CHECKING LAMBDA VALUES
//pneg[ar]=exp(-lambda*veca[ar]*4);// the prob of neg at each age group and year
//pprimsero1[ar]=(exp(-3*lambda*veca[ar]))*(1-exp(-lambda*veca[ar]));
//pprimsero2[ar]=pprimsero1[ar];
//pprimsero3[ar]=pprimsero1[ar];
//pprimsero4[ar]=pprimsero1[ar];

    exposure1=sum(lambda[1,vecj[startLambda:endLambda]]); // expectation of DV1 exposure
    exposure2=sum(lambda[2,vecj[startLambda:endLambda]]);
    exposure3=sum(lambda[3,vecj[startLambda:endLambda]]);
    exposure4=sum(lambda[4,vecj[startLambda:endLambda]]);

    pneg[ar]=exp(-exposure1-exposure2-exposure3-exposure4);//prob of neg each age group & year
    pprimsero1[ar]=(1-exp(-exposure1))*exp(-exposure2-exposure3-exposure4);
    pprimsero2[ar]=(1-exp(-exposure2))*exp(-exposure1-exposure3-exposure4);
    pprimsero3[ar]=(1-exp(-exposure3))*exp(-exposure1-exposure2-exposure4);
    pprimsero4[ar]=(1-exp(-exposure4))*exp(-exposure1-exposure2-exposure3);
        
    psec[ar]= 1-pneg[ar]-pprimsero1[ar]-pprimsero2[ar]-pprimsero3[ar]-pprimsero4[ar];// the probabilities of each secondary infection at each age group and year
        
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
  
// CHECKING LAMBDA VALUES
// lambda ~ beta(2,5);

  for (i in 1:4){
    for (j in 1:(l_vecj/step_year)){
      lambda[i,j] ~ beta(2,5); // lambda is a matrix => vector distribution
    }
  }
    
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

