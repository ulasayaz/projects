function [xhat] = PrimalDual(Phi,z,tau,lambda)

maxItersPD=30;
theta=0.9; beta=0.25; sigma=0.25; 
p=Phi.frame_length;
% start value:
xhat = zeros(p,1);      % coefficients in Phi-frame
xbar=zeros(p,1);


% start value:
xi=zeros(Phi.N-[0 0 1]); %get the size of the data cube (1 less frame)

wait_handle = waitbar(0,['PrimalDual ',...
    '(itmax=',num2str(maxItersPD),') ...']);

iter=1;
% START Primal-Dual iteration
while  (iter <= maxItersPD)  %&& (normresrel > OptTol)
    
    waitbar(iter / maxItersPD);
    
    % reconstruct x in the naturaL basis from xhat:
    x=diff(Phi.synthesize(xbar),1,3); % D o Phi^(-1) composition
    
    w=xi+sigma*x;
    
    % apply P_K*(w) soft-thresholding (ST) with parameter tau*(1-lambda)
    if isreal(w(1))
        xi  = w- (abs(w) > tau*(1-lambda)) .* ...
            (w - tau*(1-lambda).*sign(w));
    else
        xi  = w- (abs(w) > tau*(1-lambda)) .* ...
            (w - tau*(1-lambda).*angle(w)); % ST also valid for the complex case.
    end
    
    xhat0=xhat;
    wbar= xhat-beta*Phi.analyze(diff_adjoint(xi)); % Phi o D*
    
    
    % apply P_H(z)(wbar) soft-thresholding with parameters tau,lambda,beta
    wbar=(z+wbar)/2; thres=beta*tau*lambda/2;
    
    if isreal(wbar(1))
        xhat  = (abs(wbar) > thres) .* ...
            (wbar - thres.*sign(wbar));
    else
        xhat  = (abs(wbar) > thres) .* ...
            (wbar - thres.*angle(wbar)); % ST also valid for the complex case.
    end
    
    xbar=xhat+theta*(xhat-xhat0);
    
    %norm(xhat)
    
    iter=iter+1;
    
end % END Primal-Dual iteration

close(wait_handle );