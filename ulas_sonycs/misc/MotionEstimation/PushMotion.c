#include <stdio.h>
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
/* y=PushMotion(x,Dx,Dy) pushes values of matrix x to target of vector field Dx, Dy */
/* difference to PushMotion: does NOT set occluded cells (marked by nan) to 0 */
/* if required, this step must be done after this routine has finished */
/* remark: pull motion on adjoint vector field is equivalent to push motion on original*/
{
    double *Dx, *Dy, *x, *y;
    int i, j, i2, j2, M, N, idx;    
	double nan;  // machine dependent constants    

    if ((nrhs != 3) || (nlhs !=1)) {
        printf("Error! Expecting exactly 3 rhs and 1 lhs argument!\n");
        return;
    }
	 
    x  = mxGetPr(prhs[0]);
    Dx = mxGetPr(prhs[1]);
	Dy = mxGetPr(prhs[2]);	
    M  = mxGetM(prhs[0]);
    N  = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    y = mxGetPr(plhs[0]);    
    /* initialize output y to x */   
    for(i=0; i<M*N; i++) {                
		y[i]=x[i];        
    }
    
    for(i=0; i<M; i++) {        
        for(j=0; j<N; j++) {
			idx=(j*M)+i; // source cell
		    if ( (Dx[idx]!=0) || (Dy[idx]!=0) ) {
				i2=i+Dy[idx]; // target y-ccord (rows)
				j2=j+Dx[idx]; // target x-coord (columns)
				if ( (i2>=0) && (i2<M) && (j2>=0) && (j2<N) ) {
					y[(j2*M)+i2]=x[idx]; // push from source to to target cell
				}			
			}
         }   
      
    }
} /* end mexFunction */
