#include <Rmath.h>
#include "pomp.h"

//double expit (double x) {return 1.0/(1.0+exp(-x));}

typedef enum {
  SUSCEPTIBLE = 1,
  INFECTED = 2,
  CONVALESCENT = 3,
} type_t;

typedef struct{
  type_t type;
  int n_in;
  int n_out;
  int edge_out[5];
  int edge_in[5];
  int serotype;
} node_t;

typedef struct{
  int from;
  int to;
  double rate;
  int serotype;
} edge_t;

#define LOGBETA_SD     	(p[parindex[0]]) // environmental stochasticity SD in transmission rate
#define LOGMU          	(p[parindex[1]]) // mortality rate
#define LOGBETA 				(p[parindex[2]]) // transmission rate 1
#define LOGGAMMA      	(p[parindex[3]]) // recovery rate 1
#define LOGDELTA      	(p[parindex[4]]) // convalescence rate 1
#define LOGTHETA  			(p[parindex[5]]) // interaction during infection 1
#define LOGCHI  				(p[parindex[6]]) // interaction after recovery 1
#define LOGOMEGA				(p[parindex[7]]) // environmental reservoir

 // pars for observation model
#define LOGITRHO1				(p[parindex[8]]) // reporting rate of the observation model
#define LOGITRHO2				(p[parindex[9]]) // reporting rate of the observation model

#define N								(x[stateindex[80]]) // population
#define W								(x[stateindex[81]]) // white noise
#define CASES1 					(x[stateindex[82]]) // all cases1
#define CASES2 					(x[stateindex[83]]) // all cases2
#define CASES3					(x[stateindex[84]]) // all cases3
#define CASES4					(x[stateindex[85]]) // all cases4
#define CASESP1 				(x[stateindex[86]]) // primary cases1
#define CASESP2 				(x[stateindex[87]]) // primary cases2
#define CASESP3					(x[stateindex[88]]) // primary cases3
#define CASESP4					(x[stateindex[89]]) // primary cases4

#define Y1							(y[obsindex[0]])
#define Y2							(y[obsindex[1]])
#define Y3							(y[obsindex[2]])
#define Y4							(y[obsindex[3]])
#define YP1							(y[obsindex[4]])
#define YP2							(y[obsindex[5]])
#define YP3							(y[obsindex[6]])
#define YP4							(y[obsindex[7]])

 // vars
#define N_NODES         80  // total #nodes
#define S_NODES         16  // suscepitble #nodes
#define I_NODES         32  // infected #nodes
#define C_NODES         32  // convalescent #nodes

#define N_EDGES         225  // total #edges
#define B_EDGES         1 // births
#define D_EDGES         80 // deaths
#define I_EDGES         32 // infections
#define R_EDGES         32 // recoverys
#define C_EDGES         32 // convalescents
#define S_EDGES         48 // seroconversions


void poisson_rmeasure (double *y, double *x, double *p, 
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
	double lambda1, lambda2, lambda3, lambda4, lambdap1, lambdap2, lambdap3, lambdap4, rho1, rho2;
	rho1 = expit(LOGITRHO1);
	rho2 = expit(LOGITRHO2);
	lambda1 = rho1*CASES1;
	lambda2 = rho1*CASES2;
	lambda3 = rho1*CASES3;
	lambda4 = rho1*CASES4;
	lambdap1 = rho2*CASESP1;
	lambdap2 = rho2*CASESP2;
	lambdap3 = rho2*CASESP3;
	lambdap4 = rho2*CASESP4;
	GetRNGstate();

	Y1 = (R_FINITE(lambda1)) ? rpois(lambda1) : R_NaReal;
	Y2 = (R_FINITE(lambda2)) ? rpois(lambda2) : R_NaReal;
	Y3 = (R_FINITE(lambda3)) ? rpois(lambda3) : R_NaReal;
	Y4 = (R_FINITE(lambda4)) ? rpois(lambda4) : R_NaReal;
	YP1 = (R_FINITE(lambdap1)) ? rpois(lambdap1) : R_NaReal;
	YP2 = (R_FINITE(lambdap2)) ? rpois(lambdap2) : R_NaReal;
	YP3 = (R_FINITE(lambdap3)) ? rpois(lambdap3) : R_NaReal;
	YP4 = (R_FINITE(lambdap4)) ? rpois(lambdap4) : R_NaReal;

	PutRNGstate();

}

void poisson_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
	double tol = 1.0e-17;
	double lambda1, lambda2, lambda3, lambda4, rho1, rho2;
	double lik1, lik2, lik3, lik4;

	rho1 = expit(LOGITRHO1);
	rho2 = expit(LOGITRHO2);
	lambda1 = rho1*CASES1;
	lambda2 = rho1*CASES2;
	lambda3 = rho1*CASES3;
	lambda4 = rho1*CASES4;
	//Rprintf("%lg,%lg,%lg,%lg\n",lambda1,lambda2,lambda3,lambda4);
	//Rprintf("%lg,%lg,%lg,%lg\n",Y1,Y2,Y3,Y4);
	if(R_FINITE(lambda1) && R_FINITE(lambda2) && R_FINITE(lambda3) && R_FINITE(lambda4)){
		lik1 = dpois(Y1, lambda1, 0) + tol;
		lik2 = dpois(Y2, lambda2, 0) + tol;
		lik3 = dpois(Y3, lambda3, 0) + tol;
		lik4 = dpois(Y4, lambda4, 0) + tol;
		*lik = lik1*lik2*lik3*lik4;
	} else{
		*lik = tol*tol*tol*tol;
	}
	//Rprintf("%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",*lik,CASES1,Y1,lik1,CASES2,Y2,lik2);
	if(give_log) *lik = log(*lik);
	
}
void poisson_rmeasure_cum (double *y, double *x, double *p, 
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
	double lambda1, lambda2, rho1, rho2;
	rho1 = expit(LOGITRHO1);
	rho2 = expit(LOGITRHO2);
	lambda1 = rho1*CASES1;
	//lambda2 = rho2*CASES2;
	GetRNGstate();
	if(R_FINITE(lambda1)){
		Y1 = rpois(lambda1);
		Y2 = Y1;
	} else{
		Y1 = R_NaReal;
		Y2 = R_NaReal;
	}
	PutRNGstate();

}

void poisson_dmeasure_cum (double *lik, double *y, double *x, double *p, int give_log,
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
	double tol = 1.0e-17;
	double lambda1, lambda2, rho1, rho2, lambda;
	double lik1, lik2;

	rho1 = expit(LOGITRHO1);
	rho2 = expit(LOGITRHO2);
	lambda1 = rho1*CASES1;
	lambda2 = rho2*CASES2;
	lambda = lambda1 + lambda2;
	// Yall = Y1 + Y2;
//	Rprintf("%lg,%lg,%lg,%lg\n",CASES1,CASES2,lambda1,lambda2);
	if(R_FINITE(lambda1)){
		lik1 = dpois(Y1, lambda, 0) + tol;
		//lik2 = dpois(Y2, lambda2, 0) + tol;
		*lik = lik1;
	} else{
		*lik = tol;
	}
//	Rprintf("%lg,%lg,%lg,%lg,%lg\n",*lik,CASES1,Y1,lik1,Y2);
	if(give_log) *lik = log(*lik);
	
}

#undef Y1
#undef Y2
#undef Y3
#undef Y4
#undef YP1
#undef YP2
#undef YP3
#undef YP4

void msi_euler_simulate (double *x, const double *p, 
			 const int *stateindex, const int *parindex, const int *covindex,
			 int covdim, const double *covar, 
			 double t, double dt)
{
  
  double mu, beta, gamma, delta, theta, chi, omega;
  double beta_var, beta_sd, dW;
  int i, j, k;
          
  double lambda_1 = 0.0;
  double lambda_2 = 0.0;
  double lambda_3 = 0.0;
  double lambda_4 = 0.0;
  
  double xtheta_1 = 0.0;
  double xtheta_2 = 0.0;
  double xtheta_3 = 0.0;
  double xtheta_4 = 0.0;

  double rate[N_EDGES];
  double trans[N_EDGES];
  
  int nbranch[15] = {4,3,3,3,3,2,2,2,2,2,2,1,1,1,1}; // no. of main branches
  int brancht[32] = {1,2,3,4,5,6,7,5,8,9,6,8,10,7,9,10,11,12,11,13,12,13,11,14,12,14,13,14,15,15,15,15}; // where the main branches go
  int nbranchsc[28] = {3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1}; // no. of sc branches
  int branchsct[48] = {7,10,13, 3,10,13, 3,6,13, 3,6,9, 14,16, 11,16, 11,13, 15,17, 8,17, 8,14, 12,16,7,16,7,11,11,13,6,13,6,10, 13,11,12,9,10,8,9,5,7,4,5,3}; // where sc branches go
  int st[32] = {4,3,2,1, 3,2,1, 4,2,1, 4,3,1, 4,3,2, 2,1, 3,1, 3,2, 4,1, 4,2, 4,3, 1,2,3,4}; // serotypes of the main branches
  int st2[48] = {3,2,1, 4,2,1, 4,3,1, 4,3,2, 2,1, 3,1, 3,2, 2,1, 4,1, 4,2, 3,1, 4,1, 4,3, 3,2, 4,2, 4,3, 1,2,1,3,2,3,1,4,2,4,3,4}; // sertypes of the sc branches
  
  // untransform the parameters
  mu = exp(LOGMU);
  beta = exp(LOGBETA);
  gamma = exp(LOGGAMMA);
  delta = exp(LOGDELTA);
  theta = exp(LOGTHETA);
  chi = exp(LOGCHI);
  omega = exp(LOGOMEGA);
  beta_sd = exp(LOGBETA_SD);
  beta_var = beta_sd*beta_sd;
  
  // test
  if (!(R_FINITE(mu)) ||
          !(R_FINITE(beta)) ||
          !(R_FINITE(gamma)) ||
          !(R_FINITE(delta)) ||
          !(R_FINITE(theta)) ||
          !(R_FINITE(chi)) ||
          !(R_FINITE(omega)) ||
          !(R_FINITE(N)) ||
          !(R_FINITE(W)))
  return;


  if (beta_sd > 0.0){
    dW = rgamma(dt/beta_var,beta_var); // gamma noise, mean=dt, variance=(beta_sd^2 dt)
    if (!(R_FINITE(dW))) return;
  } else {			// environmental noise is OFF
    dW = dt;
  }
  
  // initialize all 80 nodes
  node_t node[N_NODES];
  for(i=0; i<N_NODES; i++){
    if(i <S_NODES){ //susceptibles
      node[i].type=SUSCEPTIBLE;
      node[i].serotype=0;
    } else if (i >=S_NODES+I_NODES) { // convalescents
      node[i].type= CONVALESCENT;
      node[i].serotype=st[i-S_NODES-I_NODES];
    } else { // infecteds
      node[i].type=INFECTED;
      node[i].serotype=st[i-S_NODES];
    }
    node[i].n_in=0;
    node[i].n_out=0;
  }	
  
  // initialize all 225 edges
  edge_t edge[N_EDGES];
  // births 
  edge[0].from=-1;
  edge[0].to=0;
  edge[0].rate=mu;
  edge[0].serotype=0;
  // deaths
  for(j=1; j<B_EDGES+D_EDGES; j++){
    edge[j].from=j-1;
    edge[j].to=-1;
    edge[j].rate=mu;
    edge[j].serotype=0;
  }
  
  int index=0;
  for(i=0; i<15; i++){
    for(j=0; j<nbranch[i]; j++){
      // infections
      edge[D_EDGES+B_EDGES+index].from=i;
      edge[D_EDGES+B_EDGES+index].to=16+index;
      edge[D_EDGES+B_EDGES+index].rate=0.0;
      edge[D_EDGES+B_EDGES+index].serotype=st[index];            
      // progression to convalescent
      edge[D_EDGES+B_EDGES+I_EDGES+index].from=16+index;
      edge[D_EDGES+B_EDGES+I_EDGES+index].to=48+index;
      edge[D_EDGES+B_EDGES+I_EDGES+index].rate=gamma;
      edge[D_EDGES+B_EDGES+I_EDGES+index].serotype=st[index]; 
      // recovery from convalescent
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+index].from=48+index;
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+index].to=brancht[index];
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+index].rate=delta;
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+index].serotype=st[index]; 
      index=index+1;
    }
  }
  // seroconversion

  int index2=0;
  for(i=0; i<28; i++){
    for(j=0; j<nbranchsc[i]; j++){
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+C_EDGES+index2].from=48+i;
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+C_EDGES+index2].to=48+i+branchsct[index2];
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+C_EDGES+index2].rate=0.0;
      edge[D_EDGES+B_EDGES+I_EDGES+R_EDGES+C_EDGES+index2].serotype=st2[index2];
      index2=index2+1;
    }
  } 

  int totinf1=0;
  int totinf2=0;
  int totinf3=0;
  int totinf4=0;
  
  for(i=S_NODES; i<S_NODES+I_NODES; i++){
    if(node[i].type==INFECTED){
      // Rprintf("i=%d, nos=%d;\n",i, x[stateindex[i]]);
      if(node[i].serotype==1){
        totinf1 = totinf1 + x[stateindex[i]];
      }else if(node[i].serotype==2){
        totinf2 = totinf2 + x[stateindex[i]];
      }else if(node[i].serotype==3){
        totinf3 = totinf3 + x[stateindex[i]];
      }else if(node[i].serotype==4){
        totinf4 = totinf4 + x[stateindex[i]];
      }
    }
  }

  // Rprintf( "T1=%d, T2=%d, T3=%d, T4=%d;\n", totinf1, totinf2, totinf3, totinf4);

  // sero-specific forces of infections
  lambda_1 = (beta*totinf1/N + omega)*dW/dt;
  lambda_2 = (beta*totinf2/N + omega)*dW/dt;
  lambda_3 = (beta*totinf3/N + omega)*dW/dt;
  lambda_4 = (beta*totinf4/N + omega)*dW/dt;

  xtheta_1 = theta*lambda_1;
  xtheta_2 = theta*lambda_2;
  xtheta_3 = theta*lambda_3;
  xtheta_4 = theta*lambda_4;

  // Rprintf( "L1=%lf, L2=%lf, L3=%lf, L4=%lf;\n", lambda_1, lambda_2, lambda_3, lambda_4);

  // Enchancement
  double mult =0.0;

  for(j=B_EDGES+D_EDGES; j<B_EDGES+D_EDGES+I_EDGES; j++){
      mult = 1.0;
    if(j>=B_EDGES+D_EDGES+4){
      mult = chi;
    }
    if(edge[j].serotype==1){
      edge[j].rate=mult*lambda_1;  
    } else if(edge[j].serotype==2){
      edge[j].rate=mult*lambda_2;  
    } else if(edge[j].serotype==3){
      edge[j].rate=mult*lambda_3;  
    } else if(edge[j].serotype==4){
      edge[j].rate=mult*lambda_4;  
    }  
     
  }	
  // Silent seroconversion
  for(j=N_EDGES-S_EDGES; j<N_EDGES; j++){
    if(edge[j].serotype==1){
      edge[j].rate=xtheta_1;  
    } else if(edge[j].serotype==2){
      edge[j].rate=xtheta_2;  
    } else if(edge[j].serotype==3){
      edge[j].rate=xtheta_3;  
    } else if(edge[j].serotype==4){
      edge[j].rate=xtheta_4;  
    }   
  }

  // node-edge graph
  for(j=0; j<N_EDGES; j++){
    if(edge[j].from >-1){
      //Rprintf("j=%d, ef=%d, st=%d, er=%lf; \n",j, edge[j].from, edge[j].serotype, edge[j].rate);
      node[edge[j].from].edge_out[node[edge[j].from].n_out] = j;
      node[edge[j].from].n_out += 1;
    }
    if(edge[j].to >-1){
      // Rprintf("j=%d, et=%d; \n",j, edge[j].to);
      node[edge[j].to].edge_in[node[edge[j].to].n_in] = j;
      node[edge[j].to].n_in += 1;
    }
  }

  /*
  for(i=0; i<80; i++){
    Rprintf("\n i=%d, eo=%d, nt=%d, nst=%d, ", i, node[i].n_out, node[i].type, node[i].serotype);
    for(j=0; j<node[i].n_out; j++){
      Rprintf("e:%d, et:%d, nst%d; ", node[i].edge_out[j], edge[node[i].edge_out[j]].to, node[edge[node[i].edge_out[j]].to].serotype);
    }
  }
   */
  rate[0] = mu*N;
  int count=1;
  for(i=0; i<N_NODES; i++){
    for(j=0; j<node[i].n_out; j++){
      rate[count] = edge[node[i].edge_out[j]].rate;
      count +=1;
    }
  }   

  trans[0] = rpois(rate[0]*dt);

  count=0;
  for(i=0; i<N_NODES; i++){
    reulermultinom(node[i].n_out, x[stateindex[i]], &rate[count+1], dt, &trans[count+1]);
    if(node[i].type==INFECTED){
      for(k=0; k<node[i].n_out; k++){
		    if(edge[node[i].edge_out[k]].serotype==1){
		      CASES1 += trans[1+count+k];
					// Rprintf("edge=%d\n",node[i].edge_out[k]);
					if(node[i].edge_out[k] < B_EDGES+D_EDGES+I_EDGES +4){
						CASESP1 += trans[1+count+k];
					}
		    }else if(edge[node[i].edge_out[k]].serotype==2){
		      CASES2 += trans[1+count+k];
					if(node[i].edge_out[k] < B_EDGES+D_EDGES+I_EDGES +4){
						CASESP2 += trans[1+count+k];
					}
		    }else if(edge[node[i].edge_out[k]].serotype==3){
		      CASES3 += trans[1+count+k];
					if(node[i].edge_out[k] < B_EDGES+D_EDGES+I_EDGES +4){
						CASESP3 += trans[1+count+k];
					}
		    }else if(edge[node[i].edge_out[k]].serotype==4){
		      CASES4 += trans[1+count+k];
					if(node[i].edge_out[k] < B_EDGES+D_EDGES+I_EDGES +4){
						CASESP4 += trans[1+count+k];
					}
		    }
			}
    }
    count = count+node[i].n_out;
  }


  // bits of code that spits out data to construct graph
	/*
  for(j=0; j<80; j++){
    Rprintf("%d, %d , %d \n", j, node[j].type, node[j].serotype);
  }
  
  for(j=0; j<225; j++){
    Rprintf("%d, %d, %d, %d \n", j, edge[j].from, edge[j].to, edge[j].serotype);
  }
      	*/
  /*
  for(j=0; j<80; j++){
    Rprintf("%d, %d , %d \n", j, node[j].type, node[j].n_out);
  }
  */	
  /*  
  for(j=0; j<225; j++){
    Rprintf("j=%d, rate=%lf, trans=%lf; \n", j, rate[j], trans[j]);
  }
  */

  // add up all the births
  N +=  trans[0];
  x[stateindex[0]] += trans[0];

  count=0;
  for(i=0; i<80; i++){
    N += -trans[count+1];
    for(j=0; j<node[i].n_out; j++){
      x[stateindex[i]] += - trans[count+1];
      if(edge[node[i].edge_out[j]].to>-1){
        x[stateindex[edge[node[i].edge_out[j]].to]] += trans[count+1];
      }
      count += 1;
    }
  }

  if (beta_sd > 0.0) {
    W += (dW-dt)/beta_sd;	// mean zero, variance = dt
  }
}

#undef LOGBETA_SD    
#undef LOGMU 
#undef LOGBETA         
#undef LOGGAMMA       
#undef LOGDELTA     
#undef LOGTHETA
#undef LOGCHI
#undef LOGITRHO1
#undef LOGITRHO2
#undef LOGOMEGA

#undef N
#undef W
#undef CASES1
#undef CASES2
#undef CASES3
#undef CASES4
#undef CASESP1
#undef CASESP2
#undef CASESP3
#undef CASESP4

#undef N_NODES
#undef S_NODES
#undef I_NODES
#undef C_NODES
#undef N_EDGES
#undef B_EDGES        
#undef D_EDGES         
#undef I_EDGES        
#undef R_EDGES         
#undef C_EDGES       
#undef S_EDGES
