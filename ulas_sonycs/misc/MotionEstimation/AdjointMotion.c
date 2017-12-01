#include <stdio.h>
#include "mex.h"
#include <math.h>  

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
/* [Ax,Ay]=AdjointMotion(Dx,Dy) adjoint motion to Dx, Dy */
{
    double *Dx, *Dy, *Ax, *Ay;
    int i, j, i2, j2, M, N, idx, idx2;    
    double Anorm, Dnorm;
    double nan, inf;  // machine dependent constants
    bool Aisnan;

    if ((nrhs != 2) || (nlhs !=2)) {
        printf("Error! Expecting exactly 2 rhs and 2 lhs arguments!\n");
        return;
    }
	 
    Dx = mxGetPr(prhs[0]);
	Dy = mxGetPr(prhs[1]);	
    M  = mxGetM(prhs[0]);
    N  = mxGetN(prhs[0]);
    
    nan=mxGetNaN();   // get machine representation for nan (not a number)
    inf=mxGetInf();
   
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    Ax = mxGetPr(plhs[0]);    
    plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    Ay = mxGetPr(plhs[1]);
    
    // points where no vectors end get nan:
    // hence initialize output matrices to nan:
    for(i=0; i<M*N; i++) {                
		Ax[i]=nan;        
        Ay[i]=nan;     
    }
     
    /* loop over matrix  */
    for(i=0; i<M; i++) {        
        for(j=0; j<N; j++) {
			idx=(j*M)+i;   // present 1d index
            i2=i+Dy[idx];  // shift in y direction (rows)
            j2=j+Dx[idx];  // shift in x direction (columns)                             
            if ( (i2>=0) && (i2<M) && (j2>=0) && (j2<N) ) {
                Dnorm=fabs(Dx[idx])+fabs(Dy[idx]);   
                idx2=(j2*M)+i2; // shifted 1d index
                Anorm=fabs(Ax[idx2])+fabs(Ay[idx2]);
                Aisnan=( (mxIsNaN(Ax[idx2])) && (mxIsNaN(Ay[idx2])) ) ;
                if ( Aisnan || (Anorm==0) ) {
                    // nan and 0 can be overwritten
                    Ax[idx2]=-Dx[idx];
                    Ay[idx2]=-Dy[idx];
                } else if (Dnorm>0) {
                    // mark crossings of 2 non-zero motions (adjoint vf is ambiguous)
                    Ax[idx2]=nan;
                    Ay[idx2]=2;
                    // assert: !Aisnan && Anorm!=0
                }
            }
 
//             Dnorm=fabs(Dx[idx])+fabs(Dy[idx]);       
//             if ( (Dnorm>Anorm) || (Aisnan) ) {
//                 // nan must always be overwritten, non-zero shift must not be overwritten by 0-shift
//                 Ax[idx2]=-Dx[idx];
//                 Ay[idx2]=-Dy[idx];
//             }
      
        }
    }
} /* end mexFunction */