L=16;
s2=Signal2D.make_fromImage('cameraman.bmp');
videoflag=2; motionunit=[1,1]; %  [0,0];  % [0.1,0.1]
signal=Signal3D.make_CyclefromLowerDimSig(s2,videoflag,L,motionunit);

Phi=Wavelet3D_mlab(signal);

B=HartleyTrafo3(signal);

c=8; % RL=60;nodes=B.ftrafo.nodesOnRadialLines(RL,c);
nodes=B.nodesDefault(c);
CS=CompressiveSensing();
CS.set_original(signal);

% CS.sim_SampledONS(Phi,B,nodes); % instead
% set system
CS.set_APhiNodes(Phi,B,nodes);
% simulate measurement vector:
CS.set_y(CS.A.sim_measurement());

%used properties: CS.y, CS.Phi, CS.A
% paramaters: DR=lambda, tau. PD= beta, sigma, theta
% maxItersDR, maxItersPD

lambda=0.5; tau=1;
maxItersDR=30;

epsilon=0;
normy=max(eps,norm(CS.y));

if ~CS.Phi.isTrafodone
    CS.Phi.dec; % to get the same size for the output Phi.dec(Phi.deepestlev)
end
p=CS.Phi.frame_length;
% start value:
xhat = zeros(p,1);      % coefficients in Phi-frame

% one needs pseudo inverse, backslash (QR) does not work!
% computational expensive for large frame.data !!!
Ainv=[];
if ~isempty(CS.A.frame)
    Ainv=pinv(CS.A.frame.data);
end

wait_handle = waitbar(0,['DouglasRachford ',...
    '(itmax=',num2str(maxItersDR),') ...']);

CS.solutionProps=struct;

iter=1;
% START Douglas-Rachford iteration
while  (iter <= maxItersDR)  %&& (normresrel > OptTol)
    
    waitbar(iter / maxItersDR);
    
    % Since xhat are the coefficients in the frame Phi, e.g. wavelet,
    % reconstruct x in the naturaL basis from xhat:
    x=CS.Phi.synthesize(xhat);
    
    % residuum is measured in the natural basis:
    if isempty(CS.A.frame)
        % fast op
        res=CS.y-CS.A.recAnalysis(x);
    else
        % matrix op
        res=CS.y-CS.A.frame.data*x(:);
    end
    normres=norm(res);
    % normresrel0=normresrel;
    normresrel=normres/normy;
    
    % projection in the natural basis:
    rproj = max(1 -epsilon/normres,0).*res;
    
    % transform residuum backwards to Phi-frame
    if isempty(Ainv)
        rprojhat=CS.A.recSynthesis(rproj);
    else
        rprojhat=Ainv* rproj(:);
    end
    % correction in Phi-frame
    corr=CS.Phi.analyze(rprojhat);
    
    % result of proximal mapping P_F of first operator F
    % xhat1=P_F(xhat)=xhat+corr:
    xhat1 = xhat + corr(:);
    
    % xhatdiff is its argument
    xhatdiff = 2*xhat1 - xhat;
    xhat = PrimalDual(CS.Phi,xhatdiff,tau,lambda)-corr;
    
    norm(xhat1)
    normresrel
    
    iter=iter+1;
end % END while

close(wait_handle );

CS.x=CS.xorig.make_like();
CS.x.replace_signal(xhat1);

normresrel=norm(res)/normy;
CS.solutionProps.numIters = iter-1;
CS.solutionProps.normresrel=normresrel;
CS.solutionProps.l1normSolution=norm(xhat1,1);
CS.solver_paramstr=['\gamma=',num2str(tau,'%3.1g'),...
    ', \epsilon=',num2str(epsilon,'%3.1g')];

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





