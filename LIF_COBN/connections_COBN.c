/* 
 * @file connections_COBN.c
 * @brief mex code to simulate a recurrent random network with excitatory 
 *        and inhibitory LIF CONDUCTANCE-BASED neurons. 
 *          The network is fully described in the paper 
 *          "Comparison of the dynamics of neural interactions between
 *          current-based and conductance-based integrate-and-fire 
 *          recurrent networks" written by S.Cavallari, S.Panzeri 
 *          and A.Mazzoni and published in Frontiers in Neural Circuits 
 *          (2014), 8:12. doi:10.3389/fncir.2014.00012.
 *          Please cite this paper if you use the code.
 * @author Stefano Cavallari
 * @date 2 2014
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
// Needed for building mex file
#include <mex.h>
#include <matrix.h>

/* Random number generator routine. To generate real random numbers 0.0-1.0
 * Should be seeded with a negative integer
*/
#include "ran1.h"  

/* gaussian variate generator routine. To generate random number along a normal distribution */
#include "gasdev.h" 


 
 /* Conventions for variable names:
 * -------------------------------
 * - prefix "e" stands for "excitatory"
 * - prefix "i" stands for "inhibitory"
 * - prefix "x" stands for "external"
 *
 * - prefix "a" stands for "AMPA"
 * - prefix "g" stands for "GABA"
 * 
 * - prefix "N"   stands for "number of..."
 * - prefix "tot" stands for "total"
 * 
 * - prefix "e2e" stands for "excitatory to excitatory"
 * - prefix "e2i" stands for "excitatory to inhibitory"
 * - prefix "x2e" stands for "external to excitatory"
 * - prefix "e2i" stands for "external to inhibitory"
 *
 * - "T"  stands for greek letter "tau"
 * 
 * - "nrn"  stands for "neurons"
 * - "FR"   stands for "firing-rate"
 * - "SP"   stands for "spike"
 * - "SPC"  stands for "spike-count"
 *
 * - "rec" stands fo recurrent
 * - "ext" stands fo external
*/


/*
 * Main entry called from Matlab
 * @param nlhs number of left-hand-side arguments
 * @param plhs pointers to mxArrays in the left-hand-side
 * @param nrhs number of right-hand-side arguments
 * @param prhs pointers to mxArrays in the right-hand-side
 *
 * This is the entry point of the function
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Variables declaration ----------------------------------------------
     
    mxArray *tmp;
    
    /* Number of neurons: */
    mwSize  eNnrn;   /* number of excitatory neurons */
    mwSize  iNnrn;   /* number of inhibitory neurons */
    mwSize  totNnrn; /* total number of neurons */
    mwIndex nrn;     /* neuron index */
    mwIndex nrn2;    /*    "     "   */
    
	double RAND_MAX_double; /* max value returned by function random*/
    long seed1;           /* first random seed (it defines connections) */
    long seed2;           /* second random seed (it defines poisson variate) */
    
    
    // Connections:
    double   p;
    mwSize  *A;             /* matrix storing neurons connections */
    mwSize  *Ncon;          /* number of outgoing connections of each neuron */
    
    mwSize *syntemp;
    mwIndex kn, ip;
    int avNcon_outgo;
    int sdNcon_outgo;
    
    /* Check of the number of inputs. For v1, net_COBN is the only input */
    if(nrhs != 1)
        mexErrMsgTxt("Not enough input arguments.");

    /* Extract the network parameters from net_COBN */

    /* Number of neurons */
    tmp   = mxGetField(prhs[0], 0, "eNnrn");
    eNnrn = *mxGetPr(tmp);
    
    tmp   = mxGetField(prhs[0], 0, "iNnrn");
    iNnrn = *mxGetPr(tmp);
    
    totNnrn = eNnrn + iNnrn;
    
    /* Connection probability */
    tmp = mxGetField(prhs[0], 0, "p");
    p   = *mxGetPr(tmp);
   
        
    /* On the use of the random seeds:
     * -------------------------------
     * The two random seeds are used in the following way:
     * - the first random seed, seed1, affects the network configuration.
     * - the second random seed, seed2, affects the generation of Poisson
     *   random variates, which determine the input to each neuron
     * When seed1==0 or seed2==0 is passed as input, time(NULL) is used instead to generate seed1 and seed2.
    */
     
    tmp   = mxGetField(prhs[0], 0, "SEED_connections" );
    seed1 = *mxGetPr(tmp);
    if(seed1==0) seed1 = time(NULL);

    tmp   = mxGetField(prhs[0], 0, "SEED_poisson" );
    seed2 = *mxGetPr(tmp);
    if(seed2==0) seed2 = time(NULL);
    
    if (seed1<0 || seed2 <0)
    {
        mexErrMsgTxt("seed1 and seed2 must be positive integer\n");
    }
    
    /* Confirm input parameters */
	mexPrintf("C code using input parameters received:\n");
    mexPrintf("                   eNnrn = %d\n", eNnrn);
    mexPrintf("                   iNnrn = %d\n", iNnrn);
    mexPrintf("                       p = %f\n", p);
    mexPrintf("seed1 (SEED_connections) = %d\n", seed1);
    mexPrintf("seed2 (SEED_poisson)     = %d\n\n", seed2);

    
    /* Allocating matrices */
    A       = mxCalloc(totNnrn * totNnrn, sizeof(mwSize)); /* array with the connections of each neuron */
    Ncon    = mxCalloc(totNnrn, sizeof(mwSize));           /* array with the number of outgoing connections of each neuron */
    
    /* Connections --------------------------------------------------------
     * To generate the random connections we make use of the seed that has
     * been provided by the user (if it has been provided), seed1, so that the same
     * network configuration can be used repeatedly.
    */
        
    syntemp = mxCalloc(totNnrn * totNnrn, sizeof(mwSize)); /* auxiliary array to build the connections of the network */
	RAND_MAX_double = (double) RAND_MAX;
    seed1 = -seed1; /* seed for ran1 must be a negative integer */

    for(nrn=0; nrn<totNnrn; nrn++) { /* nrn is the postsynaptic neuron */
        //recurrent AMPA connections entering the nrn-th neuron
        for(nrn2=0;nrn2<eNnrn;nrn2++) syntemp[nrn2]=0; /* initializing the syntemp array */
        do
        {   /* Avoid endless loop */
            ip = p*eNnrn + (int)(0.5+gasdev(&seed1)*sqrt((float)(p*eNnrn))); /* number of AMPA connections entering the nrn-th neuron */
        } while (ip > eNnrn);
        /* The mean number of AMPA connections entering each neuron is (p*eNnrn) */
        /* The variance of the number of AMPA connections entering each neuron is (p*eNnrn) */
        for(nrn2=0; nrn2<ip; nrn2++) {     
            do
            {
                kn = ran1(&seed1)*eNnrn; /* randomly selecting an excitatory presynaptic neuron */
            } while ((kn == nrn) || (kn >= totNnrn));   /* Avoid self connection & keep in range */
            if (syntemp[kn]==0) {
                syntemp[kn]=1;
                A[Ncon[kn] + kn*totNnrn] = nrn; /* connection: kn -> nrn */
                Ncon[kn]++; /* number of outgoing connections of the kn-th neuron */
            }
            else{   
                nrn2--;
            }
        }
        //GABA connections entering the nrn-th neuron
        for(nrn2=0;nrn2<iNnrn;nrn2++) syntemp[nrn2]=0; /* initializing the syntemp array */
        do
        {   /* Avoid endless loop */
            ip = p*iNnrn + (int)(0.5+gasdev(&seed1)*sqrt((float)(p*iNnrn))); /* number of GABA connections entering the nrn-th neuron */
        } while (ip > iNnrn);
        /* The mean number of GABA connections entering each neuron is (p*iNnrn) */
        /* The variance of the number of GABA connections entering each neuron is (p*iNnrn) */
        for(nrn2=0; nrn2<ip; nrn2++) {
            do
            {
                kn = eNnrn + ran1(&seed1)*iNnrn; /* randomly selecting an inhibitory presynaptic neuron */
            } while ((kn == nrn) || (kn >= totNnrn));   /* Avoid self connection & keep in range */
            if (syntemp[kn - eNnrn]==0) {
                syntemp[kn - eNnrn]=1;
                A[Ncon[kn] + kn*totNnrn] = nrn; /* connection: kn -> nrn */
                Ncon[kn]++; /* number of outgoing connections of the kn-th neuron */
            }
            else{   
                nrn2--;
            }
        }
    }
    mxFree(syntemp);


    // check of the average number of outgoing connections ---------------------
    avNcon_outgo=0;
    sdNcon_outgo=0;
    for(nrn=0; nrn<totNnrn; nrn++)
    {    
        avNcon_outgo += Ncon[nrn];
    }

    avNcon_outgo = avNcon_outgo/totNnrn;
    for(nrn=0; nrn<totNnrn; nrn++)
    {
        sdNcon_outgo += pow(Ncon[nrn]-avNcon_outgo,2);
    }
    sdNcon_outgo = sqrt(sdNcon_outgo/(totNnrn-1));
   
    mexPrintf("Av number of outgoing connections per neuron = %i \nStd of the number of outgoing connections per neuron = %i\n", avNcon_outgo, sdNcon_outgo);
    
    //--------------------------------------------------------------------------
    

    /* All done. Now copy only A to outputs (Ncon can be regenerated as required  */

	// A
    plhs[0] = mxCreateNumericMatrix(totNnrn, totNnrn, mxINT64_CLASS, mxREAL);
    mwSize *AcopyPr = (mwSize*)mxGetPr(plhs[0]);
    for(nrn=0; nrn<totNnrn*totNnrn; nrn++)
    {
        AcopyPr[nrn] = (mwSize)A[nrn];
    }

    mxFree(Ncon);
    mxFree(A);

    mexPrintf("connections_COBN: Array generation done.\n\n");  
} // main end

