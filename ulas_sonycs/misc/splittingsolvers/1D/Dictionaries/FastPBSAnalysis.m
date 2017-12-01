function c = FastPBSAnalysis(x, par1, par2, par3)
% FastPBSAnalysis -- Fast Wavelet(periodized,bi-orthogonal,symmetric) 
%			Analysis Operator
%  Usage:
%	c = FastPBSAnalysis(x, L, qmf, dqmf)
%  Inputs
%    x		1-d signal; length(x) = 2^J
%    L      	Coarsest scale (2^(-L) = scale of V_0); L << J;
%    qmf    	quadrature mirror filter
%    dqmf   	dual quadrature mirror filter (symmetric, dual of qmf)
%  Outputs
%    C      	1-d wavelet transform of x
%
%  Description
%    1. qmf filter may be obtained from MakeBSFilter   
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_PBS
%
%  See Also
%    FastPBSSynthesis, FWT_PBS, IWT_PBS, MakeBSFilter
%

n = length(x);
c = FWT_PBS(x, par1, par3, par2);
c = c(:);
