function x = tvmin(y,lambda,numdeep,lmin,lmax)
% Fast TV-limization: Chambolle-Darbon algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||_2 + lambda || grad(x) ||_1 
%
% Usage
%	x = tvmin(y,lambda,numdeep,lmin,lmax)
% Input
%	y           nxn image
%	lambda	    regularization parameter
%   	numdeep	    depth of the dyadic search: 2^numdeep levels are actually
%		    computed (hence the precision is (lmax-lmin)*2^(-numdeep))
%	lmin, lmax  the solution is truncated to lmin/lmax: for an optimal result with 
%		    no truncation they must be set to the actual min/max of the original G.
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

if ~exist('lmin'), 	lmin = min(y(:)); end
if ~exist('lmax'), 	lmax = max(y(:)); end
if ~exist('numdeep'), 	numdeep = 8; 	  end
if numdeep <=0 | numdeep > 16,
	error('Bad depth');
	return;
end

x = tvmin_mex(y,lambda,numdeep,lmin,lmax);

	

