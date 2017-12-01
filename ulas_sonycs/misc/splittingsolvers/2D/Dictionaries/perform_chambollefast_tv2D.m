function x = perform_chambollefast_tv2D(y,lambda,tvoptions)
% Fast TV-limization: Chambolle-Darbon algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1 
%
% Usage
%	x = perform_chambollefast_tv2D(y,lambda,tvoptions)
% Input
%	y           nxn image
%	lambda	    regularization parameter
%	tvoptions   a structure containing options for TV prox subiterations:
%   			numdeep	    depth of the dyadic search: 2^numdeep levels are actually
%		    	            computed (hence the precision is (lmax-lmin)*2^(-numdeep))
%			lmin, lmax  the solution is truncated to lmin/lmax: for an optimal result with 
%		                    no truncation they must be set to the actual min/max of the original G.
% Outputs
%	 x	    solution of the problem
%
% Fast TV-minimization by min graph-cuts
%  TV4 = with nearest neighbours interaction
%  TV8 = with also next nearest
%  TV16 = with 16 neighbours
%  
%  Based on the code by 
%      A. Chambolle and J. Darbon: On total variation
%      minimization and surface evolution using parametric maximum flows,
%      preprint (2008).
%  Their code implements Dorit Hochbaum's algorithm:
%     D. S. Hochbaum: An efficient algorithm for image segmentation,
%     Markov random fields and related problems. J. ACM, 48(4):686--701,
%     2001.	
%  
%

%if ~exist('lambda'),	 lambda = 1;		end
%if ~exist('tvoptions'),  
%	tvoptions.numdeep = 4;
%	tvoptions.lmin = min(y(:));
%	tvoptions.lmax = max(y(:));
%else
%	if ~isfield(tvoptions,'numdeep'),    tvoptions.numdeep = 8;	 end
%	if ~isfield(tvoptions,'lmin'),       tvoptions.lmin = min(y(:)); end
%	if ~isfield(tvoptions,'lmax'),       tvoptions.lmax = max(y(:)); end
%	if tvoptions.numdeep <=0 | tvoptions.numdeep > 16,
%		error('Bad depth');
%		return;
%	end
%end


n = floor(sqrt(prod(size(y))));
y = reshape(y, n, n);

x = tvmin_mex(y,lambda,tvoptions.numdeep,tvoptions.lmin,tvoptions.lmax);

x = x(:);

	

