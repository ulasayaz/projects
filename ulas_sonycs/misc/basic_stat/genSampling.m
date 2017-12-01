function [mask,stat,actpctg] = genSampling(pdf,params)
% [mask,stat,N] = ~(pdf[,params}) sampling pattern (with minimal interference
% if params.iter>1).
%
% a monte-carlo algorithm to generate a sampling pattern with
% minimum peak interference. The number of samples will be
% sum(pdf) *(1 +-tol)
%
%	pdf - probability density function to choose samples from
%	params.iter - number of tries
%	params.tol  - the relative deviation from the desired number of samples in samples
%
% returns:
%	mask - sampling pattern
%	stat - vector of min interferences measured each try
%	actpctg    - actual undersampling factor
%
%	(c) Michael Lustig 2007
%   G.Troll April 2014:   -- change of signature (add. parameters moved to params)
%                         -- extended to simple sampling (disregarding interference)
%                            if params.iter=1
%                         -- secure loop in Beroulli trials by itmaxBT
%                         -- changed tolerance tol to from absolute to relative value
%

    function BernoulliTrial()
        % find a mask_candidate with the required number of nodes        
        d=Inf;
        it=0;
        while d > ctol && it<=itmaxBT
            % Bernoulli trials with success prob. given by pdf
            mask0 = rand(size(pdf))<pdf;
            it=it+1;
            d0=abs(sum(mask0(:)) - K);
            if d0<d
                d=d0;
                mask_candidate=mask0;
            end
        end
    end

if nargin < 2
    params=struct;
end
if ~isfield(params,'tol') || params.tol>1
    % relative error of required sample nodes
    params.tol=0.01;
end
if ~isfield(params,'iter')     
    params.iter=10;
end

K = sum(pdf(:));  % number of sampled nodes (positive outcomes)
ctol= max(1,K*params.tol); % max. absolute deviation of number of sampled nodes

%h = waitbar(0);

pdf(pdf>1) = 1;

minIntr = Inf;
mask_candidate = false(size(pdf));
stat=zeros(1,params.iter);

itmaxBT=max(20,ceil(256^2/numel(pdf)));  % increase itmax for short signals


if params.iter>1  % minimize interference
    for n=1:params.iter
        
        % draw random mask with required number of sampled nodes
        BernoulliTrial();
        
        % determine interference pattern
        TMP = ifftn(mask_candidate./pdf);
        % TMP(1)=0; figure, imagesc(ifftshift(abs(TMP))); colorbar;
        % s=Signal2D(mask_candidate./pdf-1); f=FourierTrafos2(s); f.ifft2;f.show_resultd;
        mv=max(abs(TMP(2:end)));
        if mv < minIntr
            % choose mask with minimal peak(max) interference
            minIntr = mv;
            mask = mask_candidate;
        end
        stat(n) = mv;
        %waitbar(n/params.iter,h);
    end
else   % determine sample with required number of sampled nodes (disregarding interference)
    BernoulliTrial();
    mask=mask_candidate;
end

actpctg = sum(mask(:))/numel(mask);

%close(h);
end

