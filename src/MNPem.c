/******************************************************************
  This file is a part of MNP: R Package for Estimating the 
  Multinomial Probit Models by Kosuke Imai, Jordan R. Vance, and 
  David A. van Dyk.
  Copyright: GPL version 2 or later.
*******************************************************************/

#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"

void cMNPem(int *piNSamp, int *piNMaxit, int *piNEit, int *piNInc, int
	    *piNCov, int *piNDim, double *pdX, int *y, double *pdbeta,
	    double *pdSigma, int *verbose, double *pdStore, double
	    *pdNorm, double *pdBridge){

  /* paramerters from R */
  int n_samp = *piNSamp;   /* sample size */
  int n_gen = *piNMaxit;   /* max number of iterations */
  int n_cov = *piNCov;     /* The number of covariates */
  int n_dim = *piNDim;     /* The number of indentifiable dimension (p-1) */
  int M = *piNEit;         /* maximum number of MC draws for W */
  int inc = *piNInc;       /* increment for the number of MC draws */
  int m, m_pre;            /* number of MC draws at current iteration */
    
  double *beta, *beta_pre; /* The model parameters */
  double alpha2;           /* unidentified parameter */
  double **Sigma, **Sigma_pre; /* The unidentifiable variance matrix */
  double **X;              /* The covariates and outcome var */
  double *bridge, **Xbeta, **Xbeta_pre; /* used for bridge sampling */

  /* total number of parameters */
  int n_par = n_cov+n_dim*(n_dim+1)/2; 

  /* model parameters */
  double **W, ***Wstore;     /* The missing data! */
  double **EW, ***EWW;       /* conditional expectations */ 
  double cmean;  	     /* The conditional mean for Wij */
  double cvar;               /* The conditional variance for Wij */
  double maxw=0.0, minw=0.0; /* used for sampling W */
  double ***PerSig;          /* Permuted Sigma for sampling W */
  int kindx, lindx;          /* Indexes for the permutation of Sigma */
  double **Xstar;            /* L*X where L is the Chol factor of SigInv */
  double **SigInv;           /* The Inverse of Sigma: unidentified version */
  double **SS;	             /* Sum of squares for sweep operator */
  int i, j, k, l, main_loop, mc_loop; /* used for loops */

  /* temporay storages */
  int itemp;
  double dtemp, *vtemp, *vtemp1;
  double **mtemp,**mtemp1;

  /** get random seed **/
  GetRNGstate();

  /** defining vectors and matricies **/
  W = doubleMatrix(n_samp, n_dim);
  EW = doubleMatrix(n_samp, n_dim);
  EWW = doubleMatrix3D(n_samp, n_dim, n_dim);
  Wstore = doubleMatrix3D(M, n_samp, n_dim);
  X = doubleMatrix(n_samp*n_dim, n_cov+1);
  Xbeta = doubleMatrix(n_samp, n_dim);
  Xbeta_pre = doubleMatrix(n_samp, n_dim);
  SS = doubleMatrix(n_cov+1, n_cov+1);
  Sigma = doubleMatrix(n_dim, n_dim);
  Sigma_pre = doubleMatrix(n_dim, n_dim);
  SigInv = doubleMatrix(n_dim, n_dim);
  beta = doubleArray(n_cov);
  beta_pre = doubleArray(n_cov);
  bridge = doubleArray(n_samp);
  vtemp = doubleArray(n_dim-1);
  vtemp1 = doubleArray(n_dim);
  Xstar = doubleMatrix(n_samp*n_dim, n_cov+1);
  mtemp = doubleMatrix(n_dim, n_dim);
  mtemp1 = doubleMatrix(n_dim, n_cov);
  PerSig = doubleMatrix3D(n_dim, n_dim, n_dim);

  /** Packing Y, X, beta, Sigma  **/
  itemp = 0;
  for (k = 0; k < n_cov; k++) 
    for (j = 0; j < n_dim; j++) 
      for (i = 0; i < n_samp; i++) X[i*n_dim+j][k] = pdX[itemp++];

  itemp = 0;
  for (j=0; j<n_cov;j++) {
    beta[j]=pdbeta[itemp++];
    pdStore[j]=beta[j];
  }

  itemp = 0;
  for (k = 0; k < n_dim; k++) 
    for (j = 0; j < n_dim; j++) 
      Sigma[j][k] = pdSigma[itemp++];
  dinv(Sigma,n_dim,SigInv); 
  alpha2=Sigma[0][0];

  itemp = 0;
  for (k = 0; k < n_dim; k++) 
    for (j = 0; j < n_dim; j++) 
      if(j<=k) {
	pdStore[n_cov+itemp]=Sigma[j][k];
	itemp++;
      }

  /** initialize W **/
  for(j=0;j<n_samp;j++)     
    for(k=0;k<n_dim;k++) 
      if(y[i]==(j+1))
	W[j][k]=0.1;
      else
	W[j][k]=-0.1;

  /*** PX-MCECM algorithm! ***/
  for(main_loop=1; main_loop<=n_gen; main_loop++){
    /*** E STEP ***/
    /** permutation of Sigma **/
    for(j=0;j<n_dim;j++){
      kindx = 0;
      for(k=0;k<n_dim;k++){ 
	lindx=0;
	for(l=0;l<n_dim;l++){ 
	  if(j!=k)
	    if(j!=l)
	      PerSig[j][k+kindx][l+lindx]=Sigma[k][l];
	    else{
	      lindx=-1;
	      PerSig[j][k+kindx][n_dim-1]=Sigma[k][l];
	    }
	  else{
	    kindx=-1;
	    if(j==l){
	      lindx=-1;
	      PerSig[j][n_dim-1][n_dim-1]=Sigma[j][j];
	    }
	    else
	      PerSig[j][n_dim-1][l+lindx]=Sigma[k][l];
	  }
	}
      }
      dinv(PerSig[j],n_dim,PerSig[j]);
    }
    
    /** Truncated Multinomial Normal Sampling for W **/
    for(i=0;i<n_samp;i++) {
      for(j=0;j<n_dim;j++) { /* initialize the matrices */
	EW[i][j]=0;
	Xbeta[i][j]=0;
	for(k=0;k<n_dim;k++) EWW[i][j][k]=0;
	for(k=0;k<n_cov;k++) Xbeta[i][j]+=X[i*n_dim+j][k]*beta[k];
      }
      m_pre = m; 
      m = imin2(M, main_loop*inc); 
      for(mc_loop=0;mc_loop<m;mc_loop++) {
	for(j=0;j<n_dim;j++){
	  maxw=0.0; itemp=0;
	  for(k=0;k<n_dim;k++) {	  		  		  
	    if(j!=k) {
	      if(maxw<W[i][k]) maxw=W[i][k];
	      vtemp[itemp]=W[i][k];
	      for(l=0;l<n_cov;l++) 
		vtemp[itemp]-=X[i*n_dim+k][l]*beta[l];
	      itemp++;
	    }
	  }
	  /* conditional mean and variance */
	  cmean=Xbeta[i][j];
	  cvar=1/PerSig[j][n_dim-1][n_dim-1];
	  for(k=0;k<(n_dim-1);k++) 
	    cmean-=PerSig[j][n_dim-1][k]*vtemp[k]*cvar;
	  /* sampling each W[i][j] conditionally on the other elements */
	  if(y[i]==(j+1)) 
	    W[i][j]=TruncNorm(maxw,cmean+100*sqrt(cvar),cmean,cvar); 
	  else
	    W[i][j]=TruncNorm(cmean-100*sqrt(cvar),maxw,cmean,cvar);
	}
	/* conditional expectations */
	for(j=0;j<n_dim;j++) {
	  Wstore[mc_loop][i][j]=W[i][j];
	  EW[i][j]+=(W[i][j]*sqrt(alpha2)/m);
	  for(k=0;k<n_dim;k++)
	    EWW[i][j][k]+=(W[i][j]*W[i][k]*alpha2/m);
	}
      }
      for(j=0;j<n_dim;j++) 
	X[i*n_dim+j][n_cov]=EW[i][j];
    }
    /** bridge sampling part I **/
    if(main_loop>1) {
      pdBridge[main_loop-1]=0;
      for(i=0;i<n_samp;i++) {
	dtemp=0;
	for(j=0;j<m_pre;j++)
	  dtemp+=sqrt(exp(dMVN(Wstore[j][i], Xbeta_pre[i], Sigma_pre,
			       n_dim, 1) - dMVN(Wstore[j][i], Xbeta[i], Sigma, n_dim, 1))); 
	bridge[i]/=dtemp;
	pdBridge[main_loop-1]+=log(bridge[i]);
      }
    } 

    /*** CM STEP ***/
    /* multiply X and W by the Inverse of the Cholesky factor */
    dcholdc(SigInv,n_dim,mtemp);
    for(i=0;i<n_samp*n_dim;i++)
      for(j=0;j<=n_cov;j++) Xstar[i][j]=0;
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++) 
	  for(l=0;l<=n_cov;l++)
	    Xstar[i*n_dim+k][l]+=mtemp[j][k]*X[i*n_dim+j][l];
    
    /* construct SS matrix for SWEEP */
    for(j=0;j<=n_cov;j++)
      for(k=0;k<=n_cov;k++) SS[j][k]=0;
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++)
	for(k=0;k<=n_cov;k++)
	  for(l=0;l<=n_cov;l++) 
	    SS[k][l]+=Xstar[i*n_dim+j][k]*Xstar[i*n_dim+j][l];
    
    /* SWEEP to update beta and X*beta */
    for(j=0;j<n_cov;j++) SWP(SS,j,n_cov+1);
    for(j=0;j<n_cov;j++) {
      beta_pre[j]=beta[j];
      beta[j]=SS[j][n_cov];
    }
    for(i=0;i<n_samp;i++)
      for(j=0;j<n_dim;j++) {
	Xbeta_pre[i][j]=Xbeta[i][j];
      }

    /* update Sigma given beta */
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++) {
	Sigma_pre[j][k]=Sigma[j][k];
	Sigma[j][k]=0;
      }
    for(i=0;i<n_samp;i++) {
      for(j=0;j<n_dim;j++) 
	vtemp1[j]=0;
      for(j=0;j<n_dim;j++) {
	for(k=0;k<n_dim;k++) 
	  Sigma[j][k]+=EWW[i][j][k];
	for(k=0;k<n_cov;k++) {
	  mtemp1[j][k]=2*EW[i][j]*beta[k]; /* E(W)*beta^T */
	  vtemp1[j]+=X[i*n_dim+j][k]*beta[k]; /* X*beta */
	}
      }
      for(j=0;j<n_dim;j++) 
	for(k=0;k<n_dim;k++) 
	  Sigma[j][k]+=vtemp1[j]*vtemp1[k];
      for(j=0;j<n_cov;j++) 
	for(k=0;k<n_dim;k++) 
	  for(l=0;l<n_dim;l++) 
	    Sigma[k][l]-=mtemp1[k][j]*X[i*n_dim+l][j];
    }
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++)
	Sigma[j][k]/=(double)n_samp;

    dinv(Sigma, n_dim, SigInv); /* SigInv is unidentified version */
    alpha2=Sigma[0][0];
    for(j=0;j<n_cov;j++)
      beta[j]/=sqrt(alpha2);
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++)
	Sigma[j][k]/=alpha2;

    /** bridge sampling part II **/
    for(i=0;i<n_samp;i++){
      bridge[i]=0;
      for(j=0;j<m;j++)
	bridge[i]+=sqrt(exp(dMVN(Wstore[j][i], Xbeta[i], Sigma, n_dim, 1) - 
			dMVN(Wstore[j][i], Xbeta_pre[i], Sigma_pre, n_dim, 1)));
    } 

    /** storing the results and calculating the norm **/
    itemp=0; pdNorm[main_loop]=0;
    for(j=0;j<n_cov;j++) {
      pdNorm[main_loop]+=(beta[j]-beta_pre[j])*(beta[j]-beta_pre[j]);
      pdStore[main_loop*n_par+itemp]=beta[j];
      itemp++;
    }
    for(j=0;j<n_dim;j++)
      for(k=0;k<n_dim;k++)
	if(j<=k) {
	  pdNorm[main_loop]+=(Sigma[j][k]-Sigma_pre[j][k])*(Sigma[j][k]-Sigma_pre[j][k]);
	  pdStore[main_loop*n_par+itemp]=Sigma[j][k];
	  itemp++;
	}        
    pdNorm[main_loop]=sqrt(pdNorm[main_loop]);

    /** printing **/
    if(*verbose) {
      printf("%5d%14g%14g", main_loop, pdNorm[main_loop], pdBridge[main_loop-1]);
      for(j=0;j<n_cov;j++)
	printf("%14g", beta[j]);
      for(j=0;j<n_dim;j++)
	for(k=0;k<n_dim;k++)
	  if(j>=k)
	    printf("%14g", Sigma[j][k]);
      printf("\n");
      fflush(stdout);
    }
  } /* end of PX-MCECM */
  
  /** write out the random seed **/
  PutRNGstate();

  /** freeing memory **/
  FreeMatrix(W, n_samp);
  FreeMatrix(EW, n_samp);
  Free3DMatrix(EWW, n_samp, n_dim);
  FreeMatrix(X, n_samp*n_dim);
  FreeMatrix(SS, n_cov+1);
  FreeMatrix(Sigma, n_dim);
  FreeMatrix(SigInv, n_dim);
  free(beta);
  free(vtemp);
  free(vtemp1);
  FreeMatrix(Xstar, n_samp*n_dim);
  FreeMatrix(mtemp, n_dim);
  FreeMatrix(mtemp1, n_dim);
  Free3DMatrix(PerSig, n_dim, n_dim);
} /* main */
 
