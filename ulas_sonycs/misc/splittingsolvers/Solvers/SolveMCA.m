function part=SolveMCA(y,dict,pars1,pars2,pars3,itermax,expdecrease,lambdaStop,epsilon,verbose)
% MCA: Morphological Component Analysis of a signals and images using highly redundant dictionaries and sparstity promoting penalties.
%	   The optimization pb is solved using a modified version of the BCR algorithm.
%	   MCA solves the following optimization problem with varying lambda >= 0,
%		(part_i,for all i) = argmin Sum_i || \Phi^T_i part_i ||_p + lambda * ||y - Sum_i part_i||_2^2
%			Sparsity is promoted by the lp norm. Ideally the l_0 norm but this is relaxed to the l_1 norm.	
%				p = 1 (l_1 norm: Soft thesholding as a solution).
%				p = 0 (l_0 norm: difficult but approximated with a Hard thresholding).
%	   Each component is supposed to be sparsely described in its corresponding dictionary \Phi.
%  Usage:
%    part=SolveMCA(y,dict,pars1,pars2,pars3,itermax,expdecrease,lambdaStop,verbose)
%  Inputs:
%    y	     		input signal or image, nxn or nx1, n = 2^J
%    dict		Names of dictionaries for each part (see directory Dictionary for available transforms)
%    parsi		Parameters of dictionary (built using the MakeList function)
%    itermax		Nb of relaxation iterations
%    expdecrease	Exponential/Linear decrease of the regularization parameter
%    lambdaStop	 	Stop criterion, the algorithm stops when lambda <= lambdaStop
%    epsilon	        Stops if residual l_2 norm <= epsilon
%    verbose		
%  Outputs:
%    parti 		Estimated ith semantic component
%
%  Description
%    The dictionaries and their parameters can be built using the MakeList function.
%

% Initializations.
	if ndims(y) > 2,
		disp('Only 1D or 2D arrays are handled.');
		return;
	end
	[n,m] = size(y);
	J   = nextpow2(n);
	numberofdicts = LengthList(dict);
	if min(n,m) == 1,
		dim = 1;
		part=zeros(numberofdicts,1,n);
	else
		dim = 2;
		if n~=m,
		   disp('Only square images are handled.');
		   return;
		end
		y   = reshape(y,n,n);
		part=zeros(numberofdicts,n,n);
	end
	ind = 0;
	for nb=1:numberofdicts
	     NAME   = NthList(dict,nb);
	     PAR1   = NthList(pars1,nb);
	     PAR2   = NthList(pars2,nb);
	     PAR3   = NthList(pars3,nb);
	     DimDomain(nb) = eval(['SizeOfDict' num2str(dim) '(n, NAME, PAR1, PAR2, PAR3)']);
	     redundancy(nb) = DimDomain(nb)/prod(size(y));
	end
	thdtype='H';
	
% First pass: coeffs of the original image in each dictionary.
	coeff = eval(['FastA' num2str(dim) '(y,dict,pars1,pars2,pars3)']);

% Calculate the starting thd, which is the minimum of maximal coefficients 
% of the image in each dictionary.
	E=eval(['computeL2norms' num2str(dim) 'D(n,dict,pars1,pars2,pars3)']);
	deltamax=StartingPoint(coeff,dict,dim,pars1,pars2,pars3,DimDomain,E);
	delta=deltamax;
	if expdecrease	lambda=(deltamax/lambdaStop)^(1/(1-itermax)); % Exponential decrease.
	else	 	lambda=(deltamax-lambdaStop)/(itermax-1);     % Slope of the linear decrease. 
	end
	
% Create and return a handle on the waitbar.
	if verbose,
	   h = waitbar(0,'MCA in progress: Please wait...');
	   nbpr=ceil(sqrt(numberofdicts+2));
	end

% Start the modified Block Relaxation Algorithm.
for iter=0:itermax-1
	  % Calculate the residual.
	    residual=y-squeeze(sum(part,1));
	  
	    ind = 0;
	  % Cycle over dictionaries.
	   for nb=1:numberofdicts
	   % Update Parta assuming other parts fixed.
	   % Solve for Parta the marginal penalized minimization problem (Hard thesholding, l_1 -> Soft).
	     NAME   = NthList(dict,nb);
	     PAR1   = NthList(pars1,nb);
	     PAR2   = NthList(pars2,nb);
	     PAR3   = NthList(pars3,nb);
	     Ra=squeeze(part(nb,:,:))+residual;
	     coeffa = eval(['FastA' num2str(dim) '(Ra,NAME,PAR1,PAR2,PAR3)']);
	     if thdtype == 'H' 	coeffa  = HardThresh(coeffa,delta*E((ind+1):(ind+DimDomain(nb))));
	     else		coeffa  = SoftThresh(coeffa,delta*E((ind+1):(ind+DimDomain(nb))));
	     end
	     ind 	    = ind + DimDomain(nb);
	     part(nb,:,:)   = eval(['FastS' num2str(dim) '(coeffa,n,NAME,PAR1,PAR2,PAR3)'])/redundancy(nb);
	   end
	
	% Update the regularization parameter delta.
	    if expdecrease	delta=delta*lambda; % Exponential decrease.	
	    else		delta=delta-lambda; % Linear decrease.
	    end
	
	% Display the progress time.
	if verbose,
	    waitbar((iter+1)/itermax,h);
	    if dim==1
	      subplot(nbpr,nbpr,1);plot(y);drawnow;
	      subplot(nbpr,nbpr,2);plot(squeeze(sum(part,1)));title('\Sigma_i Part_i');drawnow;
	      for np=1:numberofdicts
	    	subplot(nbpr,nbpr,np+2);plot(squeeze(part(np,:,:)));title(sprintf('Part_%d',np));drawnow;
	      end
	    else
	      subplot(nbpr,nbpr,1);imagesc(y);axis image;rmaxis;drawnow;
	      subplot(nbpr,nbpr,2);imagesc(squeeze(sum(part,1)));axis image;rmaxis;title('\Sigma_i Part_i');drawnow;
	      for np=1:numberofdicts
	    	subplot(nbpr,nbpr,np+2);imagesc(squeeze(part(np,:,:)));axis image;rmaxis;title(sprintf('Part_%d',np));drawnow;
	      end
	    end
	end
	
	if norm(residual(:)) <= epsilon,
	  if verbose,
	  	close(h);
	  end
	  return;
	end
	% 
end

if verbose,
% Close the waitbar window
	close(h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = StartingPoint(C,dict,dim,pars1,pars2,pars3,DimDomain,E)

nbdicts = LengthList(dict);
ind = 0;
for i = 1:nbdicts,
      NAME = NthList(dict, i);
      PAR1 = NthList(pars1, i);
      PAR2 = NthList(pars2, i);
      PAR3 = NthList(pars3, i);
      tmp = C((ind+1):(ind+DimDomain(i)))./E((ind+1):(ind+DimDomain(i)));
      buf(i)=max(abs(tmp(:)));
      ind = ind + DimDomain(i);
end

%
buf=flipud(sort(buf(:)));
delta=buf(2);



