function [mv,output] = EstimateMotionVector( signal1, signal2, do_subpix, lp_norm)
% mv=~(signal1,signal2,[do_subpix]) global shift between signals
% (mv(1)=rows, mv(2)=cols);
% if do_subpix, try to reach subpixel precision (in lp-norm or cross.correlation);
% the opional return value output contains information about the
% subpixel optimization result;
% The algorithm uses phase correlation and/or lp-distance
% [according to the correlation theorem: F(Corr(f,g))=F^*(f).*F(g) ]
%

    function c=ncrosscorr(s)
       % computes the negative cross correlation at shift s
       
       % first, shift images on top of each other
       % so that the integer part of the new motion
       % vector should be zero:
       signal2b=Signal2D(Signal2D.SubPixelShifter(signal2.xn, s));
       % now compute a contour image (here diagonal signal of the
       % 1-level undecimated wavelet transform:
       w2=Wavelet2U_mlab(signal2b);
       w2.dec;
       signal2b=Signal2D(w2.detcoef(1,'d'));
       % we can compute cross-correlation only at mv2=[0,0] directly
       % instead of computing the full phase-correlation by signal1b.motionVector. 
       c=sum(signal1b.xn(:).*signal2b.xn(:));
%        % compute the phase correlation and check new motion vector:
%        [mv2,c]=signal1b.motionVector(signal2b);
%        if ~isequal(mv2,[0,0])
%            % if the integer part of the new motion
%            % vector is not zero, we have left the relevant neighbourhood
%            c=min(c,0.5*c0);  % correlation value should be smaller than the start value
%        end
       c=-c;  % go to negative for fminsearch
    end

    function c=dist(s)
       % computes the lp-distance at shift s
       
       % first, shift images on top of each other
       x2=Signal2D.SubPixelShifter(signal2.xn, s);
       % compute lp-distance of images in the canonical basis:
       c=norm(signal1.xn(:)-x2(:),lp_norm);
    end

%% main
requireOnly(isa(signal1,'Signal2D') && isa(signal2,'Signal2D'),'local',...
        'both arguments of type Signal2D');
if ~exist('do_subpix','var') || isempty(do_subpix)
    do_subpix=false;
end
if ~exist('lp_norm','var') 
    lp_norm=[];
end

% uses contour image instead of originals:
% here digonal signal of the undecimated wavelet transform:
w1=Wavelet2U_mlab(signal1);
w1.dec(1);
signal1b=Signal2D(w1.detcoef(1,'d'));

w2=Wavelet2U_mlab(signal2);
w2.dec(1);
signal2b=Signal2D(w2.detcoef(1,'d'));
[mv,corr_intpix]=signal1b.motionVector(signal2b); % corr_intpix= phase-correlation at mv
%c0=-ncrosscorr(mv);  % this c0 is cross-correlation at mv

mv0=mv;
if do_subpix
    % try to reach subpixel precision    
    mv2=[0,0];
    if isempty(lp_norm)
        [mv1,corr_subpix,exitflag,output] = fminsearch(@ncrosscorr,mv0); 
        output.corr_subpix=-corr_subpix;
    else
        [mv1,dist_subpix,exitflag,output] = fminsearch(@dist,mv0); 
        output.lp_norm=lp_norm;
        output.dist_subpix=dist_subpix;
    end
    if exitflag==1 && isequal(mv2,[0,0])
        mv=mv1;
    end
    
else
    output=[];
    mv1=[];
end
output.corr_intpix=corr_intpix;


end

