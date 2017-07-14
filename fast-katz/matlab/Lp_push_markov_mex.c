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
 * :2010-02-23: Initial coding based on katz_push_mex.c
 */

#include <mex.h>

#include <math.h>
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
        if (fabs(d[heapj2]) > fabs(d[heapj])) {
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
            if (fabs(d[right]) > fabs(d[left])) { 
                i=i+1; /* pick the smaller right child */
            }
        }
        if (fabs(d[heapk]) > fabs(d[T[i]])) {
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

const int laplacian_markov_push_history_size = 5;

/** Solve a system with the Laplacian using a Markov chain-like algorithm
 * For the system (D-A)x = ei, we can solve
 * (I - P^T)y = ei instead, where x = D^{-1} y
 * A push style algorithm handles this system assuming a small bit
 * of acceleration in the Richardson method.
 *
 * @param N, ai, aj - the parameters for the graph
 * @param ei - the non-zero right hand size index
 * @param omega - the damping parameter for sor
 * @param x the solution vector (must be allocated to size N double array)
 * @param r the residual vector (must be allocated to size N double array)
 */
void laplacian_markov_push(mwSize N, const mwIndex* ai, const mwIndex* aj, 
               mwIndex ei, double omega, 
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
    
    double st = sf_time();
    
    memset(x, 0, sizeof(double)*N);
    memset(r, 0, sizeof(double)*N);
    
    /* add ei to r, then heap it */
    r[ei] = 1;
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
        x[v] = x[v] + omega*val;
        r[v] = val - omega*val;
        val *= omega/((double)(ai[v+1]-ai[v])); /* value to propagate to neighbors */
        
        if (fabs(r[v]) > tol) {
            /* add v to the heap */
            T[n] = v;
            L[v] = n;
            n++;
            heap_up(L[v], n, T, L, r);
        }
        
        for (nzi=ai[v]; nzi<ai[v+1]; ++nzi) {
            
            /* update the residual */
            mwIndex w = aj[nzi], k;
            r[w] += val;
            
            if (fabs(r[w]) > tol && L[w]>N) {
                /* add w to the heap */
                T[n] = w;
                L[w] = n;
                n++;
                k = L[w];
                heap_up(k, n, T, L, r);
            } else if (L[w]<=N) {
                /* adjust w's position in the heap */
                k = L[w];
                k = heap_down(k, n, T, L, r);
                heap_up(k, n, T, L, r);
            }
            
            nedges += 1;
        }
        
        if (n>0) { resid = r[T[0]]; } else { resid = 0.; }
        
        if (history) {
            history[laplacian_markov_push_history_size*steps + 0] = (double)nedges;
            history[laplacian_markov_push_history_size*steps + 1] = fabs(resid);
            history[laplacian_markov_push_history_size*steps + 2] = v;
            history[laplacian_markov_push_history_size*steps + 3] = sf_time() - st;
        }
        
        if (fabs(resid)>(double)N*(1./tol)) {
            mexWarnMsgIdAndTxt("Lp_push_markov_mex:residualTooLarge",
                "after %i steps, the biggest residual is %g -- halting iteration",
                steps, resid);
            break;
        }
        
        steps++;
        if (fabs(resid) < tol) { 
            break; 
        }   
    }
    
    /* normalize */
    {mwIndex i; for (i=0; i<N; i++) { x[i] = x[i]/((double)(ai[i+1]-ai[i]));}}
    
    if (out_nedges) { *out_nedges = nedges; }
    if (out_nsteps) { *out_nsteps = steps; }
    if (out_nverts) { *out_nverts = nverts; }
}               

/* Lp_push_markov_mex
 *   [x,r,hist,nverts] = Lp_push_markov_mex(A,i,omega,tol,maxsteps)
 */
void mexFunction(int nargout, mxArray *pargout[],
                 int nargin, const mxArray* pargin[])
{
    const mxArray *argA, *argomega, *argi, *argtol, *argmaxsteps;
    mxArray *argx, *argr, *arghist, *argverts;
    double omega, tol;
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
    argi = pargin[1];
    argomega = pargin[2];
    argtol = pargin[3];
    argmaxsteps = pargin[4];
    

    /* validate arguments */
    if (!mxIsSparse(argA)) {
        mexErrMsgIdAndTxt("Lp_push_markov_mex:sparseRequired",
            "The input A argument must be a sparse matrix.");
    }
    if (!mxIsScalar(argomega)) {
        mexErrMsgIdAndTxt("Lp_push_markov_mex:scalarRequired",
            "The input alpha argument must be a scalar.");
    }
    if (!mxIsScalar(argi)) {
        mexErrMsgIdAndTxt("Lp_push_markov_mex:scalarRequired",
            "The input i argument must be a scalar.");
    }
    if (!mxIsScalar(argtol)) {
        mexErrMsgIdAndTxt("Lp_push_markov_mex:scalarRequired",
            "The input tol argument must be a scalar.");
    }
    if (!mxIsScalar(argmaxsteps)) {
        mexErrMsgIdAndTxt("Lp_push_markov_mex:scalarRequired",
            "The input argmaxsteps argument must be a scalar.");
    }
    
    /* Load the sparse matrix */
    ai = mxGetJc(argA);
    aj = mxGetIr(argA);
    n = mxGetN(argA);
    
    omega = mxGetScalar(argomega);
    i = (mwIndex)mxGetScalar(argi);
    i = i-1;
    tol = mxGetScalar(argtol);
    maxsteps = (mwSize)mxGetScalar(argmaxsteps);
    
    argx = mxCreateDoubleMatrix(n,1,mxREAL);
    argr = mxCreateDoubleMatrix(n,1,mxREAL);
    arghist = mxCreateDoubleMatrix(laplacian_markov_push_history_size,maxsteps,mxREAL);
    argverts = mxCreateDoubleMatrix(1,1,mxREAL);
    
    laplacian_markov_push(n, ai, aj, i, omega, tol, maxsteps, mxGetPr(arghist), 
        mxGetPr(argx), mxGetPr(argr), &nedges, &nverts, &nsteps);
    
    *(mxGetPr(argverts)) = (double)nverts;
    
    /* resize hist */
    mxSetN(arghist,nsteps);
    mxSetPr(arghist,mxRealloc(mxGetPr(arghist),sizeof(double)*laplacian_markov_push_history_size*nsteps));
    
    nargout = 4;
    pargout[0] = argx;
    pargout[1] = argr;
    pargout[2] = arghist;
    pargout[3] = argverts;
}
