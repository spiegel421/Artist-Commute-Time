/** @file katz_push_mex.c
 * A mex-file implementation of a push algorithm to compute Katz's scores.
 */

/*
 * David F. Gleich
 * Copyright, University of British Columbia, 2010
 */

/** 
 * History
 * -------
 * :2010-01-29: Initial coding based on katz_push.m
 */

#include <mex.h>
#include <string.h>

#include <sys/types.h>
#include <sys/timeb.h>
double sf_time()
{
#if defined(_WIN32) || defined(_WIN64)
  struct __timeb64 t; _ftime64(&t);
  return (t.time*1.0 + t.millitm/1000.0);
#else
  struct timeval t; gettimeofday(&t, 0);
  return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
#endif
}


/** Move an element up the heap until it hits the correct place.
 * @param j the index to move
 * @param n the size of the heap
 * @param T the heap items
 * @param L the location of the heap items
 * @param d the values of the heap items
 */
mwIndex heap_up(mwIndex j, mwSize n, mwIndex* T, mwIndex *L, double *d) {
    while (1) {
        if (j==0) { break; } /* the element is at the top */
        mwIndex heapj = T[j];
        mwIndex j2 = (j-1)/2;
        mwIndex heapj2 = T[j2];
        if (d[heapj2] > d[heapj]) {
            break; /* the parent is smaller, so stop */
        } else {
            /* the parent is larger, so swap */
            T[j2] = heapj; L[heapj] = j2;
            T[j] = heapj2; L[heapj2] = j;
            j = j2;
        }
    }
    return j;
}

/** Move an element down the heap until it hits the correct place.
 * @param j the index to move
 * @param n the size of the heap
 * @param T the heap items
 * @param L the location of the heap items
 * @param d the values of the heap items
 */
mwIndex heap_down(mwIndex k, mwSize n, mwIndex* T, mwIndex *L, double *d) {
    mwIndex heapk = T[k];
    while (1) {
        mwIndex i=2*(k+1)-1;
        if (i>=n) { break; } /* end of heap */
        if (i<n-1) { 
            /* pick between children (unneeded if i==heap_size-1) */
            mwIndex left=T[i];
            mwIndex right=T[i+1];
            if (d[right] > d[left]) { 
                i=i+1; /* pick the smaller right child */
            }
        }
        if (d[heapk] > d[T[i]]) {
            /* k is smaller than both children, so end */
            break;
        } else {
            T[k] = T[i]; L[T[i]]=k;
            T[i] = heapk; L[heapk] = i;
            k=i;
        }
    }
    return k;
}

void heap_print(mwSize n, mwIndex* T, mwIndex *L, double *d) {
    mwIndex i;
    for (i=0; i<n; i++) {
        mexPrintf("%6i ", T[i]);
    }
    mexPrintf("\n");
    for (i=0; i<n; i++) {
        mexPrintf("%6i ", L[T[i]]);
    }
    mexPrintf("\n");
    for (i=0; i<n; i++) {
        mexPrintf("%6.4f ", d[T[i]]);
    }
    mexPrintf("\n");
        
}

void katz_push(mwSize N, const mwIndex* ai, const mwIndex* aj, 
               double alpha, mwIndex ei,
               double tol, mwSize maxsteps, 
               double *history,
               double* x, double* r, 
               mwSize *out_nedges, mwSize *out_nverts, mwSize *out_nsteps)
{
    double val;
    /* visited array flags visited vertices */
    int *visited = mxCalloc(N,sizeof(int));
    /* create heap variables */
    mwSize n = 0;
    mwIndex *T = mxMalloc(sizeof(mwIndex)*N),
            *L = mxMalloc(sizeof(mwIndex)*N);
    /* initialize the elements of L to null-values */
    { mwIndex i; for (i=0; i<N; i++) { L[i] = N+1; } }

    mwSize nedges = 0;
    mwSize steps = 0;
    mwSize nverts = 0;
    
    
    double sumresid = 0;
    double st = sf_time();
    
    memset(x, 0, sizeof(double)*N);
    memset(r, 0, sizeof(double)*N);
    
    /* add ei to r, then heap it */
    r[ei] = 1;
    sumresid = 1.;
    T[n] = ei; 
    L[ei] = n;
    n++;
    heap_up(n-1, n, T, L, r);
    
    while (steps < maxsteps) {
        double resid;
        /* pop the top element from the heap */
        mwIndex nzi, v = T[0];
        T[0] = T[n-1];
        L[T[0]] = 0;
        L[v] = N+1; /* the null location, the order here is important! */
        n--;
        heap_down(0, n, T, L, r);

        
        if (!visited[v]) {
            visited[v] = 1;
            nverts++;
        }
        
        val = r[v];
        x[v] = x[v] + val;
        sumresid -= r[v];
        r[v] = 0.;
        val *= alpha; /* value to propagate to neighbors */
        
        for (nzi=ai[v]; nzi<ai[v+1]; ++nzi) {
            mwIndex w = aj[nzi], k;
            r[w] += val;
            sumresid += val;
            
            /* adjust w's position in the heap */
            if (L[w]>N) {
                /* add w to the heap */
                T[n] = w;
                L[w] = n;
                n++;
            }
            k = L[w];
            k = heap_down(k, n, T, L, r);
            heap_up(k, n, T, L, r);
       
            nedges += 1;
        }
        
        if (n>0) { resid = r[T[0]]; } else { resid = 0.; }
        
        if (history) {
            history[4*steps + 0] = (double)nedges;
            history[4*steps + 1] = resid;
            history[4*steps + 2] = sumresid;
            history[4*steps + 3] = sf_time() - st;
        }
        
        if (sumresid>(double)N*(1./tol)) {
            mexWarnMsgIdAndTxt("katz_push:residualTooLarge",
                "after %i steps, the residual sum is %g -- halting iteration",
                steps, sumresid);
            break;
        }
        
        steps++;
        if (resid < tol) { 
            break; 
        }   
    }
    
    if (out_nedges) { *out_nedges = nedges; }
    if (out_nsteps) { *out_nsteps = steps; }
    if (out_nverts) { *out_nverts = nverts; }
}               

/* katz_push_mex:
 *   [x,r,hist,nverts] = katz_push_mex(A,alpha,i,tol,maxsteps)
 */
void mexFunction(int nargout, mxArray *pargout[],
                 int nargin, const mxArray* pargin[])
{
    const mxArray *argA, *argalpha, *argi, *argtol, *argmaxsteps;
    mxArray *argx, *argr, *arghist, *argverts;
    double alpha, tol;
    mwIndex i;
    mwIndex *ai, *aj;
    mwSize n, maxsteps, nedges, nverts, nsteps;
    
    if (nargin < 5) {
        mexErrMsgIdAndTxt("MATLAB:nargchk:notEnoughInputs",
            "Not enough input arguments.");
    } else if (nargin > 5) {
        mexErrMsgIdAndTxt("MATLAB:nargchk:tooManyInputs",
            "Too many input arguments.");
    }
    
    argA = pargin[0];
    argalpha = pargin[1];
    argi = pargin[2];
    argtol = pargin[3];
    argmaxsteps = pargin[4];
    

    /* validate arguments */
    if (!mxIsSparse(argA)) {
        mexErrMsgIdAndTxt("katz_push_mex:sparseRequired",
            "The input A argument must be a sparse matrix.");
    }
    if (!mxIsScalar(argalpha)) {
        mexErrMsgIdAndTxt("katz_push_mex:scalarRequired",
            "The input alpha argument must be a scalar.");
    }
    if (!mxIsScalar(argi)) {
        mexErrMsgIdAndTxt("katz_push_mex:scalarRequired",
            "The input i argument must be a scalar.");
    }
    if (!mxIsScalar(argtol)) {
        mexErrMsgIdAndTxt("katz_push_mex:scalarRequired",
            "The input tol argument must be a scalar.");
    }
    if (!mxIsScalar(argmaxsteps)) {
        mexErrMsgIdAndTxt("katz_push_mex:scalarRequired",
            "The input argmaxsteps argument must be a scalar.");
    }
    
    /* Load the sparse matrix */
    ai = mxGetJc(argA);
    aj = mxGetIr(argA);
    n = mxGetN(argA);
    
    alpha = mxGetScalar(argalpha);
    i = (mwIndex)mxGetScalar(argi);
    i = i-1;
    tol = mxGetScalar(argtol);
    maxsteps = (mwSize)mxGetScalar(argmaxsteps);
    
    argx = mxCreateDoubleMatrix(n,1,mxREAL);
    argr = mxCreateDoubleMatrix(n,1,mxREAL);
    arghist = mxCreateDoubleMatrix(4,maxsteps,mxREAL);
    argverts = mxCreateDoubleMatrix(1,1,mxREAL);
    
    katz_push(n, ai, aj, alpha, i, tol, maxsteps, mxGetPr(arghist), 
        mxGetPr(argx), mxGetPr(argr), &nedges, &nverts, &nsteps);
    
    *(mxGetPr(argverts)) = (double)nverts;
    
    /* resize hist */
    mxSetN(arghist,nsteps);
    mxSetPr(arghist,mxRealloc(mxGetPr(arghist),sizeof(double)*4*nsteps));
    
    nargout = 4;
    pargout[0] = argx;
    pargout[1] = argr;
    pargout[2] = arghist;
    pargout[3] = argverts;
}
