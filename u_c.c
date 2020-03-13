/*
 * The calling syntax is:
 *
 *		u = u_c(s,t,sig,r,u0,K)
 *
 * This is a MEX-file for MATLAB.
 * 2020 Mikhailov D.E.
 *
 *========================================================*/

#include "mex.h"
#include <math.h>

/* The computational routine */
void u_c(double *u, mwSize M, mwSize N, double *s, double *sig, double *r, double *u0, double K, double hs, double ht)
{
    mwSize i;
    mwSize j;
    double *A, *B, *C, *F, *alpha, *beta;
    mxArray *pA, *pB, *pC, *pF, *palpha, *pbeta;
    pA = mxCreateDoubleMatrix(1,N,mxREAL);
    pB = mxCreateDoubleMatrix(1,N,mxREAL);
    pC = mxCreateDoubleMatrix(1,N,mxREAL);
    pF = mxCreateDoubleMatrix(1,N,mxREAL);
    palpha = mxCreateDoubleMatrix(1,N-1,mxREAL);
    pbeta = mxCreateDoubleMatrix(1,N-1,mxREAL);
    A = mxGetPr(pA);
    B = mxGetPr(pB);
    C = mxGetPr(pC);
    F = mxGetPr(pF);
    alpha = mxGetPr(palpha);
    beta = mxGetPr(pbeta);
    
    /*¬ычисление коэффициентов A,B,C */
      
    for (j=0; j<N; j++){
        if (j==0 || j==N-1){
            A[j] = 0.0;
            B[j] = 1.0;
            C[j] = 0.0;
        } else {
            A[j] = pow(s[j]*sig[j]/hs,2)/2;
            B[j] = -1/ht-(s[j]*r[j])/hs-r[j]-pow(s[j]*sig[j]/hs,2);
            C[j] = (s[j]*r[j])/hs + pow(s[j]*sig[j]/hs,2)/2;
        }
    }
    /*for (j=0; j<N; j++){
        mexPrintf("A %f", A[j]);
        mexPrintf(", B %f", B[j]);
        mexPrintf(", C %f\n", C[j]);        
    }*/
    
    /*начальное условие*/
    for (j=0; j<N; j++){
        u[j*M] = u0[j];
    }
    F[0] = K;
    F[N-1] = 0;
    alpha[0] = -C[0]/B[0];
    beta[0] = F[0]/B[0];
    for (i=0; i<M-1; i++){
        for (j=1; j<N-1; j++){
            F[j] = -u[i+j*M]/ht;
            alpha[j] = -C[j]/(A[j]*alpha[j-1]+B[j]);
            beta[j] = (F[j]-A[j]*beta[j-1])/(A[j]*alpha[j-1]+B[j]);
        }
        u[i+1+(N-1)*M] = 0;
        for (j=N-2; j>=0; j--){
            u[i+1+j*M] = alpha[j]*u[i+1+(j+1)*M]+beta[j];            
        }
    }
    /*граничные услови€
    for (i=1; i<M; i++){
        u[i] = K;
        u[i+(N-1)*M] = 0;
    } */
    
    /*for (i=0; i<M; i++) {
        for (j=0; j<N; j++){
            u[i+j*M] = j*M;
        }
    }*/
    mxDestroyArray(pA);
    mxDestroyArray(pB);
    mxDestroyArray(pC);
    mxDestroyArray(pF);
    mxDestroyArray(palpha);
    mxDestroyArray(pbeta);    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *s;  /* вектор s - цена акции */
    double *t;  /* вектор t - врем€ */
    double *sig;  /* вектор sigma - волатильность */
    double *r;  /* вектор r - ставка */
    double *u0;  /* вектор u_0 - начальное условие */
    double K;   /* K - strike */
    double *u;  /* матрица U - выход*/
            
    mwSize N;   /*длина вектора s */
    mwSize M;   /*длина вектора t */
    double hs;   /*h_s*/
    double ht;   /*h_t*/

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","6 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* получение аргументов функции */
    s = mxGetPr(prhs[0]);
    t = mxGetPr(prhs[1]);
    sig = mxGetPr(prhs[2]);
    r = mxGetPr(prhs[3]);
    u0 = mxGetPr(prhs[4]);
    K = mxGetScalar(prhs[5]);
    
    /* получение вспомогательных переменных */
    
    N = mxGetN(prhs[0]);
    M = mxGetN(prhs[1]);
    hs = s[1];
    ht = t[1];
    

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);

    /* get a pointer to the real data in the output matrix */
    u = mxGetPr(plhs[0]);

    /* call the computational routine */
    u_c(u,M,N,s,sig,r,u0,K,hs,ht);
}
