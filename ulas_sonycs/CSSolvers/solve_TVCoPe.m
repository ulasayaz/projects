%% solve_TVCoPe : 
%{
Solves a primal dual algorithm by Combettes & Pesquet (CoPe) (see final report)
min_x F(x) + G_1(L_1 x) + G_2(L_2 x) 
--(original x must be a sparse signal!)

Below we solve: 
min_x \lambda ||x||_1 + \lambda_2 || D_t \Phi^{-1} x ||_1 + 1/2 ||A Phi^{-1} x-y||_2^2
%}

%% ----- create 3d signal:

L=16;
%s2=Signal2D.make_fromImage('cameraman.bmp');
%videoflag=2; motionunit=[1.2,1.4]; %  [0,0];  % [0.1,0.1]
%signal=Signal3D.make_CyclefromLowerDimSig(s2,videoflag,L,motionunit);
L=3; fn='tennis.avi'; %fn='riksch1.avi';
signal=Signal3D.make_fromVideo(fn,L);
j=1;
s2=Signal2D(signal.xn(:,:,j)); s2.signalname=signal.signalname;
s2.colormap_active='gray';

% special case: stationary video
% signal=signal.frame(ceil(L/2)).make_VideofromSingleImage(L);


%% choose encoder and sparsifying transform

% -- case 1: 3d transform and encoder
Phi=Wavelet3D_mlab(signal);       
        
B=HartleyTrafo3(signal);

c=8; % RL=60;nodes=B.ftrafo.nodesOnRadialLines(RL,c);

% -- case2: frame-wise 2d transform and encoder
params=struct;
Phi=MultiResTransform3_Partial(signal,Wavelet2D_mlab()); 
%Phi=MultiResTransform3_Partial(signal,Curvelet2_clab());

B=FrameTrafo3_Partial(signal,HartleyTrafo2());

% B=FrameTrafo3_Partial(signal,Wavelet2D_mlab());
% B.ftrafo.set_deepestlev(1); 
% params.w1=0.995;

% -- compression rate fixed:
% c=8;

% -- fixed c, but for each frame its own nodes:
% improves reconstruction quality compared to using the same nodes
% for each frame!
c=8*ones(1,L);% RL=60;nodes=B.ftrafo.nodesOnRadialLines(RL,c);

%%

nodes=B.nodesDefault(c,params);
CS=CompressiveSensing();
CS.set_original(signal);

% CS.sim_SampledONS(Phi,B,nodes); % instead of this:
% set system

CS.set_APhiNodes(Phi,B,nodes);
% simulate measurement vector:
CS.set_y(CS.A.sim_measurement());
disp('finished');

%% CoPe reconstruction algorithm starts

%used properties: CS.y, CS.Phi, CS.A
% paramaters: gamma, lambda, lambda2, maxItersCP

% -- change lambda, lambda1, gamma for convergence
lambda=0.1; %analysis weight
lambda2=0.1; %analysis weight
% lambda2=lambda=1;
gamma=0.5;
maxItersCP=100

normy=max(eps,norm(CS.y));

if ~CS.Phi.isTrafodone
    CS.Phi.dec; % to get the same size for the output Phi.dec(Phi.deepestlev)
end

p=CS.Phi.frame_length;
% start value:
%xhat = zeros(p,1);      % coefficients in Phi-frame
xhat=CS.Phi.analyze(CS.A.recSynthesis(CS.y));
xhat=xhat(:);

v1 = zeros(p,1);
v2 = zeros(CS.Phi.N-[0 0 1]); %get the size of the data cube (1 less frame)

% one needs pseudo inverse, backslash (QR) does not work!
% computational expensive for large frame.data !!!
% Ainv=[];
% if ~isempty(CS.A.frame)
%     Ainv=pinv(CS.A.frame.data);
% end

CS.solutionProps=struct;

iter=1;

label=['CombettesPesquet (itmax=',num2str(maxItersCP),') ...'];
            multiWaitbar( label, 0); 
% START CoPe iteration
while  (iter <= maxItersCP)  %&& (normresrel > OptTol)
    
    multiWaitbar( label, iter / maxItersCP); 
    
    % Notice xhat are the coefficients in the frame Phi, e.g. wavelet.
    % calculate grad_f(xhat)
    x=CS.Phi.synthesize(xhat);
    
    % residuum is measured in the natural basis:
    if isempty(CS.A.frame)
        % fast op
        res=CS.A.recAnalysis(x)-CS.y;
        xdiff = CS.A.recSynthesis(res);
    else
        % matrix op
        res=CS.A.frame.data*x(:)-CS.y;
        xdiff = CS.A.frame.data'*res(:);
    end
    
    % transform residuum backwards to Phi-frame	where it is sparse
    xhatdiff = CS.Phi.analyze(xdiff);
    %end calculate
    % use of xhatdiff(:) is important below, in the case of MultiResTransform3_Partial
    corr=CS.Phi.analyze(diff_adjoint(v2));
    yhat = xhat - gamma*(xhatdiff(:) + v1 + corr(:));
    
    % i=1
    z1 = v1 + gamma*xhat;
    % compute proximal mapping P_G of 1-norm operator with lambda
    z1prox  = z1-(abs(z1) > lambda) .*(z1 - lambda.*sign(z1));
    v1 = v1 - z1 + z1prox + gamma*yhat;
    
    % i=2
    z2 = v2 + gamma*diff(CS.Phi.synthesize(xhat),1,3);
    % compute proximal mapping P_G of 1-norm operator with lambda2
    z2prox  = z2-(abs(z2) > lambda2) .*(z2 - lambda2.*sign(z2));
    v2 = v2 - z2 + z2prox + gamma*diff(CS.Phi.synthesize(yhat),1,3);
    % end cases
    
    % calculate grad_f(yhat)
    ybar=CS.Phi.synthesize(yhat);
    
    % residuum is measured in the natural basis:
    if isempty(CS.A.frame)
        % fast op
        resy=CS.A.recAnalysis(ybar)-CS.y;
        ydiff = CS.A.recSynthesis(resy);
    else
        % matrix op
        resy=CS.A.frame.data*ybar(:)-CS.y;
        ydiff = CS.A.frame.data'*resy(:);
    end
    
    % transform residuum backwards to Phi-frame	where it is sparse
    yhatdiff = CS.Phi.analyze(ydiff);
    %end calculate
    
    corr=CS.Phi.analyze(diff_adjoint(z2prox));
    xhat=xhat - gamma*(yhatdiff(:) + z1prox + corr(:));
    
    norm(xhat)
    %normresrel
    
    iter=iter+1;
end % END while

multiWaitbar( label, 'Close' );

CS.x=CS.xorig.make_like();
CS.x.replace_signal(xhat);

normresrel=norm(res)/normy;
CS.solutionProps.numIters = iter-1;
CS.solutionProps.normresrel=normresrel;
CS.solutionProps.l1normSolution=norm(xhat,1);
CS.solver_paramstr=['\gamma=',num2str(gamma,'%3.1g'),...
    ', \lambda=',num2str(lambda,'%3.1g'),...
    ', \lambda_2=',num2str(lambda2,'%3.1g')];

% transform from Phi-frame to natural basis:
CS.postprocess();

% compute reconstruction quality
if ~CS.xorig.isemptydata
    CS.compute_reconQuality();
end
% show result:
CS.show_recon;
% SampleTrafo22.synthesize and .recSynthesis have changed B.ts:
B.ts=CS.xorig;

%% test quality bound
cs=40;

Phi.dec(); %wavelet decomposition

test3=Phi.sparseApprox(cs); %compression and decomposition

PSNR3d=signal.PSNR(test3);
SSIM3d=signal.SSIM(test3);
mean(PSNR3d)
mean(SSIM3d)







