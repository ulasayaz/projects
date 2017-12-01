function [ x, cs ] = CS_MotionCompRecon( signal, vf, params )
% video reconstruction using motion compensation

    function disp_cond(str)
       if params.verbose
          disp(str); 
       end        
    end

if ~exist('params','var')
    params=struct;
end
if ~isfield(params,'fig') || isempty(params.fig)
    params.fig=1;
end
if ~isfield(params,'zref') || isempty(params.zref)
    params.zref=ceil(signal.size(3)/2);
end
if ~isfield(params,'verbose') || isempty(params.verbose)
    params.verbose=true;
end

disp_cond('running video reconstruction using motion compensation ...');

% best results for Phi= TV_IC_MRA3 and variable compression rates:
signal.set_zref(params.zref);
Phi= TV_IC_MRA3(signal);   % improves when c varies over frames !!!

Phi.motionCF=vf;

Phi.use_motion=true;
Phi.set_algo;
Phi.zref=params.zref;
disp_cond(['use_motion=',num2str(Phi.use_motion),', ', class(vf)]);


% -- Curvelets give further improvement BUT take longer
% Phi.set_transform(Curvelet2_clab());
% Phi.anweight1=1-1/signal.size(3);

% -- fourier ENCODER:
B=FrameTrafo3_Partial(signal,HartleyTrafo2());
L=signal.size(3);
c=8*ones(1,L);

params0=params;
params0.fig=0;
[nodes,c,cmean]=B.nodesDefault(c,params0);

% check, if nodes are ok
assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');

% create compressive-sensing object:
cs=CompressiveSensing();
cs.set_original(signal);
cs.set_methodActive(4);  % DR:1, FISTA: 2
disp_cond(cs.get_methodActive);

if L>=9
    % reduce parameters critical for convergence
    cs.params.PD.sigma=cs.params.PD.sigma/2;
    cs.params.PD.tau=cs.params.PD.tau/2;
    cs.params.PD.maxIters=500;
end

cs.fig=params.fig;
% setting epsilon >0 ---> improves rescon quality for single image but NOT for multiple frames
% cs.params.epsilon= 0.01*sqrt(signal.numel);

% simulate measurements with compression c and reconstruct:
cs.sim_SampledONS(Phi,B,nodes);
% cs.show_recon;
x=cs.x;

end

