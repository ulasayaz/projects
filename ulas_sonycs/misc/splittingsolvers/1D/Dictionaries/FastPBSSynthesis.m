function x = FastPBSSynthesis(c, n, par1, par2, par3)
% FastPBSSynthesis -- Fast Wavelet(periodized,bi-orthogonal,symmetric) 
%			Synthesis Operator
%  Usage:
%	x = FastPBSSynthesis(c, L, qmf, dqmf)
%  Inputs
%    c		1-d wavelet transform: length(wc)= 2^J.
%    L      	Coarsest scale (2^(-L) = scale of V_0); L << J;
%    qmf    	quadrature mirror filter
%    dqmf   	dual quadrature mirror filter (symmetric, dual of qmf)
%  Outputs
%    x      	1-d signal reconstructed from c
%
%  Description
%    1. qmf filter may be obtained from MakeBSFilter   
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_PBS
%
%  See Also
%    FastPBSAnalysis, FWT_PBS, IWT_PBS, MakeBSFilter
%

x = IWT_PBS(c, par1, par2, par3);
x = x(:);
