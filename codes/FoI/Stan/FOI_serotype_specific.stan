
//write function so we can sum things as sum is not built into stan


data {

// HCMC infections
int<lower=0> a;// number of age groups

int veca[a]; //vector age constitute of specific age rather than ages in order.

int datap[6,a]; ////for 2013 this will be sixgroups (neg,primary(1,2,3,4),and secondary  for all age groups)

}

parameters {

vector<lower=0, upper=1>[4] lambda;// 4 DENV serotype


}

transformed parameters {
//the probabilities of being neg, primary and secondary at each age group 
////// HCMC
// need to specify what type of data and how these are...

real pneg[a];// the probabilities of neg at each age group and year
real pprimsero1[a] ; // the probabilities of primary infection of each serotype at each age group and year
real pprimsero2[a] ;
real pprimsero3[a] ;
real pprimsero4[a] ;
real psec[a];// the probabilities of each secondary infection at each age group and year

// calculating all these probabilities
for ( ar in 1:a){

pneg[ar]=exp(-sum(lambda)*veca[ar]);// the prob of neg at each age group and year


pprimsero1[ar]=(exp(-(lambda[2]+lambda[3]+lambda[4])*veca[ar]))*(1-exp(-lambda[1]*veca[ar]));
pprimsero2[ar]=(exp(-(lambda[1]+lambda[3]+lambda[4])*veca[ar]))*(1-exp(-lambda[2]*veca[ar]));
pprimsero3[ar]=(exp(-(lambda[1]+lambda[2]+lambda[4])*veca[ar]))*(1-exp(-lambda[3]*veca[ar]));
pprimsero4[ar]=(exp(-(lambda[1]+lambda[2]+lambda[3])*veca[ar]))*(1-exp(-lambda[4]*veca[ar]));

psec[ar]= 1-pneg[ar]-pprimsero1[ar]-pprimsero2[ar]-pprimsero3[ar]-pprimsero4[ar];// the probabilities of each secondary infection at each age group and year

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


