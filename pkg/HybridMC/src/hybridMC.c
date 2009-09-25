/* *************************************************************/
/*  Hybrid Monte Carlo, with multipoint methods implemented    */
/*  C and R code by Richard D. Morey (richarddmorey@gmail.com) */
/*  For details, see Liu (2001), "Monte Carlo Strategies in    */
/*  Scientific Computing."                                     */
/*                                                             */
/***************************************************************/


#include           <R.h>
#include           <Rmath.h>  
#include           <Rdefines.h>
#include           <stdio.h>
#include           <math.h>
#include           <stdlib.h>
#include 	   <R_ext/Utils.h>



void R_CheckUserInterrupt(void);

SEXP hybridMC(SEXP ystart, SEXP nsamp, SEXP Ex, SEXP dEx, SEXP epsLo, SEXP epsRng, SEXP LFsteps, SEXP compWeights, SEXP progress, SEXP pBar, SEXP rho);
double evalEx(SEXP x, SEXP Ex, SEXP rho);
void evaldEx(SEXP x, SEXP *g, PROTECT_INDEX gpidx, SEXP dEx, SEXP rho);
SEXP MPhybridMC(SEXP ystart, SEXP nsamp, SEXP Ex, SEXP dEx, SEXP epsLo, SEXP epsRng, SEXP LFsteps, SEXP compWeights, SEXP MPwidth, SEXP MPweights, SEXP progress, SEXP pBar, SEXP rho);


/* function HybridMC arguments */

/* ystart     : initial value                                                          */
/* nsamp      : number of samples to take                                              */
/* Ex          : the "energy" function                                                  */
/* dEx         : the derivative of the "energy" function with respect to each parameter */
/* epsLo       : the lower bound of the leapfrog step size                              */
/* epsRng      : the width of the leapfrog step size range                              */
/* LFsteps     : the number of leapfrog steps to take                                   */
/* compWeights : the "weights" or "masses" of the individual variables                  */
/* progress    : How often to update the progress bar                                   */
/* pBar        : The R function to update the progress bar                              */
/* rho         : the R environment in which to execute Ex and dEx                       */ 

SEXP hybridMC(SEXP ystart, SEXP nsamp, SEXP Ex, SEXP dEx, SEXP epsLo, SEXP epsRng, SEXP LFsteps, SEXP compWeights, SEXP progress, SEXP pBar, SEXP rho)
{
	
	int j=0,i=0,N=0,S=0,s=0;
	double e=0,epsilon=0,b=0;
	double H=0,cande=0,candH=0;
	
	PROTECT_INDEX gpidx;
	
	GetRNGstate();
		
	N=length(ystart);
	S=INTEGER_VALUE(nsamp);
	
	double m[N],sumM2=0;//mStart[N];
	SEXP x, g, oldx, returnVec;
	SEXP R_fcall,sampCounter;

	PROTECT(sampCounter = NEW_INTEGER(1));
	PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(x = NEW_NUMERIC(N));
	PROTECT_WITH_INDEX(g = NEW_NUMERIC(N), &gpidx);
	PROTECT(returnVec = NEW_NUMERIC(N*S));
	PROTECT(oldx = NEW_NUMERIC(N));
	

	// make a copy of the starting value
	Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(ystart),N);

	for(s=0;s<S;s++){  // iterate through samples

	  // This really should all be done in a separate function.

	  // Copy the starting value of x for storage
	  Memcpy(NUMERIC_POINTER(oldx),NUMERIC_POINTER(x),N);
	
	  //Check to see if we need to update the progress bar
	  if(INTEGER_VALUE(progress) && !((s+1)%INTEGER_VALUE(progress))){
	    INTEGER_POINTER(sampCounter)[0]=s+1;
	    SETCADR(R_fcall, sampCounter);
	    eval(R_fcall, rho); //Update the progress bar
	  }
	  

	  // Sample new momenta, and compute the sum of squares
	  sumM2=0;
	  for(i=0;i<N;i++){
	    m[i] = rnorm(0,pow(NUMERIC_POINTER(compWeights)[i],.5));
	    sumM2+=pow(m[i],2)/(2*NUMERIC_POINTER(compWeights)[i]);
	  }
	  
	  // Copy the starting momenta for later
	  //memcpy(mStart,m,N*sizeof(double));
	  
	  // Evaluate the log-density and the derivative of the log density
	  e = evalEx(x,Ex,rho);
	  evaldEx(x, &g, gpidx,dEx, rho);
	  
	  // Compute the hamiltonian
	  H = sumM2 + e;
	  sumM2 = 0;
	  
	  
	  // Determine the value of epsilon, the time discretization parameter
	  if(epsRng==0){
	    epsilon = NUMERIC_VALUE(epsLo);
	  }else{
	    epsilon = runif(NUMERIC_VALUE(epsLo),NUMERIC_VALUE(epsLo)+NUMERIC_VALUE(epsRng));
	  }
	  
	  
	  // Leapfrog steps
	  for(i=0;i<INTEGER_VALUE(LFsteps);i++){
	    for(j=0;j<N;j++){
	      m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
	      NUMERIC_POINTER(x)[j] = NUMERIC_POINTER(x)[j] + epsilon * m[j]/NUMERIC_POINTER(compWeights)[j];
	    }
	    evaldEx(x, &g, gpidx, dEx, rho);
	    for(j=0;j<N;j++){
	      m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
	      if(i == (INTEGER_VALUE(LFsteps)-1)) sumM2+=pow(m[j],2)/(2*NUMERIC_POINTER(compWeights)[j]);
	    }
	  }
	  

	  // Determine hamiltonian for candidate
	  cande = evalEx(x,Ex,rho);
	  candH = cande + sumM2;
	  
	  

	  // Metropolis-Hastings decision
	  b=-log(runif(0,1));
	  if(!(b>(candH-H))){
	    Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(oldx),N);
	  }
	  
	  // Copy the sth sample to the vector that will be returned
	  Memcpy(NUMERIC_POINTER(returnVec)+s*N,NUMERIC_POINTER(x),N);
	  
	}
	PutRNGstate();
	UNPROTECT(6);
	
	return(returnVec);
	
}



/* function MPHybridMC arguments */
/* The MPHybridMC function implements the multipoint method of Liu (2001)


/* ystart     : initial value                                                          */
/* nsamp      : number of samples to take                                              */
/* Ex          : the "energy" function                                                  */
/* dEx         : the derivative of the "energy" function with respect to each parameter */
/* epsLo       : the lower bound of the leapfrog step size                              */
/* epsRng      : the width of the leapfrog step size range                              */
/* LFsteps     : the number of leapfrog steps to take                                   */
/* compWeights : the "weights" or "masses" of the individual variables                  */
/* progress    : How often to update the progress bar                                   */
/* pBar        : The R function to update the progress bar                              */
/* rho         : the R environment in which to execute Ex and dEx                       */ 
/* MPwidth: the width of the multipoint window                                          */
/* MPweights: the weights to be used within the multipoint window                       */

SEXP MPhybridMC(SEXP ystart, SEXP nsamp, SEXP Ex, SEXP dEx, SEXP epsLo, SEXP epsRng, SEXP LFsteps, SEXP compWeights, SEXP MPwidth, SEXP MPweights, SEXP progress, SEXP pBar, SEXP rho)
{
	
	int j=0,i=0,N=0,S=0,s=0,width=0,k=0;
	double e=0,epsilon=0,b=0;
	double H=0,candH=0;
	
	PROTECT_INDEX gpidx;
	
	GetRNGstate();
		
	N=length(ystart);
	S=INTEGER_VALUE(nsamp);
	width=INTEGER_VALUE(MPwidth);
	
	double m[N],mStart[N],sumM2=0;
	double forwardWin[width][N];
	double forwardH[width],maxForwardH=0,pCumulForward[width],pCumulBackward[width];
	double backwardH1[width],maxBackwardH2=0,maxBackwardH1=0,maxAll=0;

	double backwardH2[width];
	
	SEXP x, g, oldx, returnVec, candX;
	SEXP R_fcall, sampCounter;


	PROTECT(sampCounter = NEW_INTEGER(1));
	PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(x = NEW_NUMERIC(N));
	PROTECT_WITH_INDEX(g = NEW_NUMERIC(N), &gpidx);
	PROTECT(returnVec = NEW_NUMERIC(N*S));
	PROTECT(oldx = NEW_NUMERIC(N));
	PROTECT(candX = NEW_NUMERIC(N));

	// Copy the starting values
	Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(ystart),N);

    
	// Begin sampling S values. The value chosen at the end of the loop is used as the new starting value.
	for(s=0;s<S;s++)
	  {
	
	    Memcpy(NUMERIC_POINTER(oldx),NUMERIC_POINTER(x),N);
		
	    
	    // Check to see if we need to update the progress bar
	    if(INTEGER_VALUE(progress) && !((s+1)%INTEGER_VALUE(progress))){
	      INTEGER_POINTER(sampCounter)[0]=s+1;
	      SETCADR(R_fcall, sampCounter);
	      eval(R_fcall, rho); // update the progress bar
	    }
	    
	    // Draw the momenta and compute the sum of squares
	    sumM2=0;
	    for(i=0;i<N;i++){
	      m[i] = rnorm(0,pow(NUMERIC_POINTER(compWeights)[i],.5));
	      sumM2+=pow(m[i],2)/(2*NUMERIC_POINTER(compWeights)[i]);
	    }
	    
	    // Copy the starting momenta for later
	    memcpy(mStart,m,N*sizeof(double));
	    
	    // Compute the potential "energy" at the starting values 
	    e = evalEx(x,Ex,rho);
	    evaldEx(x, &g, gpidx, dEx, rho);
	    
	    // store the starting hamiltonian for later
	    backwardH1[0]=sumM2 + e;
	    
	    sumM2 = 0;
	    
	    // determine epsilon, the time discretization parameter
	    if(epsRng==0){
	      epsilon = NUMERIC_VALUE(epsLo);
	    }else{
	      epsilon = runif(NUMERIC_VALUE(epsLo),NUMERIC_VALUE(epsLo)+NUMERIC_VALUE(epsRng));
	    }

	    // Set the initial max H for the forward steps to the minimum possible value.
	    // We'll use the final max H as a normalizing constant
	    maxForwardH = -DBL_MAX;
	    
	    // Begin LFsteps FORWARD leapfrog steps
	    for(i=0;i<INTEGER_VALUE(LFsteps);i++)
	      {
	      sumM2=0;
	      // Move each x value
	      for(j=0;j<N;j++){
		m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
		NUMERIC_POINTER(x)[j] = NUMERIC_POINTER(x)[j] + epsilon * m[j]/NUMERIC_POINTER(compWeights)[j];
		// If they are supposed to be in the multipoint window, store them
		if(i>(INTEGER_VALUE(LFsteps) - width -1)){
		  forwardWin[i-INTEGER_VALUE(LFsteps)+width][j]=NUMERIC_POINTER(x)[j];
		}
	      }
	      
	      evaldEx(x, &g, gpidx, dEx, rho);
	      
	      for(j=0;j<N;j++){
		m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
		sumM2+=pow(m[j],2)/(2*NUMERIC_POINTER(compWeights)[j]);
	      }
	      
	      // If they are in the multipoint window, store them
	      // and check if they might be a part of the backwards window
	      // if so, store them
	      if(i>(INTEGER_VALUE(LFsteps) - width - 1)){
		forwardH[i-INTEGER_VALUE(LFsteps)+width] = evalEx(x,Ex,rho) + sumM2;
		if(forwardH[i-INTEGER_VALUE(LFsteps)+width]>maxForwardH) 
		  maxForwardH = forwardH[i-INTEGER_VALUE(LFsteps)+width];
	      }
	      if(i<(width-1)){
		backwardH1[i+1] = evalEx(x,Ex,rho) + sumM2;
	      }
	    }
	    
	    // We need to choose a value from the window.
	    pCumulForward[0] = NUMERIC_POINTER(MPweights)[0]*exp(forwardH[0]-maxForwardH);
	    for(i=1;i<width;i++){
	      pCumulForward[i] = pCumulForward[i-1] + NUMERIC_POINTER(MPweights)[i]*exp(forwardH[i]-maxForwardH);
	    }

	    
	    // Choose value randomly
	    b = runif(0,pCumulForward[width-1]);	    
	    for(i=0;i<width;i++){
	      if(b<pCumulForward[i]){
		k=i;
		break;
	      }
	    }
	    
	    // We have our candidate
	    Memcpy(NUMERIC_POINTER(candX),forwardWin[k],N);
	    

	    // We have chosen the values forwardWin[k][*]  with H = fowardH[k]
	    // We will now do the backward steps

	    // backwardH1 are the H values we've already computed
	    // they are "forward" from the starting values
	    maxBackwardH1 = -DBL_MAX;
	    for(i=0;i<(k+1);i++){
	      if(backwardH1[i]>maxBackwardH1) maxBackwardH1=backwardH1[i];
	    }
	    
	    if(k < (width-1)){ // if we are here, then we need to do a few backward steps
	      
	      maxBackwardH2 = -DBL_MAX;
	      
	      
	      //reset momenta to negative what they started at
	      sumM2=0;
	      for(i=0;i<N;i++){ 
		m[i]=-mStart[i];
		sumM2+=pow(m[i],2);
	      }
	      
	      //reset x back to starting x
	      Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(oldx),N);
	      
	      
	      // compute initial hamiltonian 
	      e = evalEx(x,Ex,rho);
	      evaldEx(x, &g, gpidx, dEx, rho);
	      
	      
	      sumM2 = 0;
	      
	      // Do our backward leapfrog steps
	      for(i=0;i<(width-k-1);i++){
		sumM2=0;
		for(j=0;j<N;j++){
		  m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
		  NUMERIC_POINTER(x)[j] = NUMERIC_POINTER(x)[j] + epsilon * m[j]/NUMERIC_POINTER(compWeights)[j];
		}
		evaldEx(x, &g, gpidx, dEx, rho);
		for(j=0;j<N;j++){
		  m[j] = m[j] - epsilon * NUMERIC_POINTER(g)[j]/2;
		  sumM2+=pow(m[j],2)/(2*NUMERIC_POINTER(compWeights)[j]);
		}
		backwardH2[width-k-2-i] = evalEx(x,Ex,rho);
		if(backwardH2[width-k-2-i]>maxBackwardH2) maxBackwardH2 = backwardH2[width-k-2-i];
	      }
	      
	      //Now we need to compute the acceptance probability
	      if(maxBackwardH2<maxBackwardH1){
		maxAll=maxBackwardH1;
	      }else{
		maxAll=maxBackwardH2;
	      }
	      
	      if(maxAll<maxForwardH) maxAll=maxForwardH;
	      
	      //maxAll is the normalizing constant we'll use
	      
	      //Renormalize all the values
	      pCumulForward[0] = NUMERIC_POINTER(MPweights)[0]*exp(forwardH[0]-maxAll);
	      for(i=1;i<width;i++){
		pCumulForward[i] = pCumulForward[i-1] + NUMERIC_POINTER(MPweights)[i]*exp(forwardH[i]-maxAll);
	      }
	      
	      pCumulBackward[0] = NUMERIC_POINTER(MPweights)[width-1]*exp(backwardH2[0]-maxAll);
	      for(i=1;i<(width-k-1);i++){
		pCumulBackward[i] = pCumulBackward[i - 1] + NUMERIC_POINTER(MPweights)[width-i-1]*exp(backwardH2[i]-maxAll);
	      }
	      
	      for(i=(width-k-1);i<width;i++){
		pCumulBackward[i] = pCumulBackward[i - 1] + NUMERIC_POINTER(MPweights)[width-i-1]*exp(backwardH1[i-width+k+1] - maxAll);
	      }
	      
	    }else{
	      //if we are here, we don't need to do any more backward steps. We have them all saved.
	      if(maxForwardH<maxBackwardH1){
		maxAll=maxBackwardH1;
	      }else{
		maxAll=maxForwardH;
	      }
	      
	      pCumulForward[0] = NUMERIC_POINTER(MPweights)[0]*exp(forwardH[0]-maxAll);
	      for(i=1;i<width;i++){
		pCumulForward[i] = pCumulForward[i-1] + NUMERIC_POINTER(MPweights)[i]*exp(forwardH[i]-maxAll);
	      }
	      
	      pCumulBackward[0] = NUMERIC_POINTER(MPweights)[width-1]*exp(backwardH1[0]-maxAll);
	      for(i=1;i<width;i++){
		pCumulBackward[i] = pCumulBackward[i-1] + NUMERIC_POINTER(MPweights)[width-i-1];
	      }
	      
	      
	    }
	    
	    
	    candH     = pCumulForward[width-1];  // These are the summed forward values
	    H         = pCumulBackward[width-1]; // These are the summed backward values
	    
	    // Metropolis-Hastings step

	    b=runif(0,1);
	    if(b>exp(candH-H)){
	      Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(oldx),N);
	    }else{
	      Memcpy(NUMERIC_POINTER(x),NUMERIC_POINTER(candX),N);
	    }
	    
	    Memcpy(NUMERIC_POINTER(returnVec)+s*N,NUMERIC_POINTER(x),N);
	    
	  }
	PutRNGstate();
	
	UNPROTECT(7); 
	
	return(returnVec);

}




/* function evalEx: evaluate the energy function */

double evalEx(SEXP x, SEXP Ex, SEXP rho)
{

	double y;
	SEXP R_fcall;

	PROTECT(R_fcall = lang2(Ex, R_NilValue));
	SETCADR(R_fcall, x);
	y = REAL(eval(R_fcall, rho))[0];
	UNPROTECT(1);

	return y;

}


void evaldEx(SEXP x, SEXP *g, PROTECT_INDEX gpidx, SEXP dEx, SEXP rho)
{

	SEXP R_fcall;

	PROTECT(R_fcall = lang2(dEx, R_NilValue));
	SETCADR(R_fcall, x);
	REPROTECT(*g = eval(R_fcall, rho), gpidx);
	UNPROTECT(1);	
	
}

