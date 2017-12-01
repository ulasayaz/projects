% checking the norm of D_bar^{-1} operator. 
dim=10;

y=eye(dim);

for i=1:dim
    for j=1:i
        y(i,j)=1;
    end
end

norm(y)

%%

%{ 
sf= Signal_Factory();
    signalBig2=Signal2D.make_fromImage('parrots.tiff');
    L=8; pansize=256*[1,1]; origin=[700,300]; motionunit=[3,5];
    signal3= sf.make_PanfromSignal2D(signalBig2,L,pansize,origin,motionunit);
    --- estimate motion
    vf=MotionVF.make_fromSignal3(signal3);
    vf.fig=5;
    vf.show_field_TV(); % all vector fields: 
    vf.show_field_TV(1); % or only the first: 
    vf.show_norm; % show norm of all vector fields: 
    vf.show_motionCompensation();
    vf.test_motionCompensation(1);
    vf.test_motionCompensation();
    %}

%% checking the linear operation for circshift

a=[1:9]';

A=reshape(a,3,3);

S=zeros(9,9);
S(1,9)=1; S(2,7)=1; S(3,8)=1; 
S(4,3)=1; S(5,1)=1; S(6,2)=1;
S(7,6)=1; S(8,4)=1; S(9,5)=1;

%% test TV_IC_MRA3_U with motion vector

fn='cameraman.bmp';
%fn='cam512.png';
L=10;
unit=1;
videoflag=2;

signal2=Signal2D.make_fromImage(fn);
%signal2.resize(128);
signal3= Signal_Factory.make_CyclefromSignal2D(signal2,videoflag,L,unit);

%{ 
y=normalize_frame(signal3.xn);
playData(y);
%}

params=struct;
% best results for Phi= TV_IC_MRA3 and variable compression rates:
Phi2= TV_IC_MRA3_U(signal3);   % improves when c varies over frames !!!
%Phi= MRA_o_TV3(signal3);   % deteriorates when c varies over frames


% -- fourier ENCODER:
B2=FrameTrafo3_Partial(signal3,HartleyTrafo2());

c=8*ones(1,L);

%RL=60;nodes=B.nodesOnRadialLines(RL);
% params.frameNumber=L;  % define for each frame new mask
% OBSERVATION: quality gets worse for varying mask
[nodes,c,cmean]=B2.nodesDefault(c,params);

% check, if nodes are ok
assert(xor(iscell(nodes.data),nodes.fL==L),'#frames is set in nodes');

% create compressive-sensing object:
CS=CompressiveSensing();
CS.set_original(signal3);
CS.set_methodActive(4);  % DR:1, FISTA: 2
disp(CS.get_methodActive);
% setting epsilon >0 ---> improves rescon quality for single image but NOT for multiple frames
% CS.params.epsilon= 0.01*sqrt(signal.numel);

% simulate measurements with compression c and reconstruct:
% parameter setting
% CS.params.epsilon=0.1;
% PD
CS.params.PD.sigma=0.05; CS.params.PD.tau=0.05; 
CS.params.PD.maxIters=500;
% DR
CS.params.DR.gamma=0.1; CS.params.DR.maxIters=200;

% SNR=[];
% SNR=3;

%decide to use motion
Phi2.motionV=unit;
Phi2.use_motion=true;

CS.sim_SampledONS(Phi2,B2,nodes);

y=normalize_frame(CS.x.xn);
playData(y);

%% test quality bound with shifted/nonshifted cameraman

fn='cameraman.bmp';
%fn='cam512.png';
L=10;
unit=0;
videoflag=2;

c=40;

signal2=Signal2D.make_fromImage(fn);
%signal2.resize(128);
signal3= Signal_Factory.make_CyclefromSignal2D(signal2,videoflag,L,unit);

Phi2= TV_IC_MRA3_U(signal3); 

Phi2.dec;

test=Phi2.sparseApprox(c);
test.graph_signal

signal3.PSNR(test)

%% wiener filtering

figure;

x = imread('cameraman.bmp');
%x = rgb2gray(x);

J = imnoise(x,'gaussian',0,0.025);
imshow(J)

K = wiener2(J,[3 3]);
figure, imshow(K)
















