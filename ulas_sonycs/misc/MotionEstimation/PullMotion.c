#include <stdio.h>
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        /* y=PullMotion(x,Ax,Ay,Q) pulls values of matrix x from target of vector field Ax, Ay */
        /* difference to PushMotion: treats occluded cells using the weight factors Q*/
        /* REMARK:
        /* --- computation of the adjoint field detects occlusions and marks them with nan;*/
        /* --- weight factors 0 will set all occluded points to 0; */
        /* --- of Q is not the 0-matrix, its support should be restricted to the occluded points; */
        /* --- pull motion on adjoint vector field is equivalent to push motion on original;*/
{
    double *Ax, *Ay, *x, *y, *Q;
    int i, j, i2, j2, M, N, idx;
    double nan;  // machine dependent constants
    bool Aisnan;
    
    if ((nrhs != 4) || (nlhs !=1)) {
        printf("Error! Expecting exactly 4 rhs and 1 lhs argument!\n");
        return;
    }
    
    x  = mxGetPr(prhs[0]);
    Ax = mxGetPr(prhs[1]);
    Ay = mxGetPr(prhs[2]);
    Q  = mxGetPr(prhs[3]);
   
    M  = mxGetM(prhs[0]);
    N  = mxGetN(prhs[0]);
    nan=mxGetNaN();   // get machine representation for nan (not a number)
    
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    // here y is initialized to 0:
    y = mxGetPr(plhs[0]);
    
    for(i=0; i<M; i++) {
        for(j=0; j<N; j++) {
            idx=(j*M)+i;  // source cell
            y[idx]=x[idx]*Q[idx];  // weight factor Q~=0 for occluded points
            /* nan marks occlusions */
            if ( !(mxIsNaN(Ax[idx])) ) {
                // move pixel by pull operation
                i2=i+Ay[idx];  // target y-coord (rows)
                j2=j+Ax[idx];  // target x-coord (columns)
                y[idx]=x[idx]*Q[idx];  // weight factor Q~=0 for occluded points
                if ( (i2>=0) && (i2<M) && (j2>=0) && (j2<N) ) {
                    y[idx]=y[idx]+x[(j2*M)+i2];  // pull back from target cell to source cell
                }
            }
        }
    }
} /* end mexFunction */
