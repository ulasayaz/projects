% Testing the 2d-Fourier transform as sparsifying transform in
% compressive sensing
% improve frequency sparsity windowing original.

if ~exist('c','var')
    c=2;
end
if ~exist('fig','var')
    fig=2;
end
if ~exist('wvn','var')
    wvn='db2';
end

signalorig=Signal2D.make_exp_symm();
%signalorig=Signal2D.make_fromImage('cameraman.bmp');
    
signal=signalorig.clone;
% set window to improve frequency sparsity:
signal.set_winfun(@chebwin);
signal=signal.apply_winfun;

Phi=HartleyTrafo2(signal);

% EITHER random sampling in standard basis (bad result) .....
B=IdTransform(signal);
% ---- OR sampling in ONS, here Walsh-Hadamard
% for better result:
% B=WalshHadamard1(signal);

nodes=B.nodesDefault(c);

cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
cs.keep_reconsparse=isa(Phi, 'FrameTrafoPair');

cs.sim_SampledONS(Phi,B,nodes);

% remove window from reconstructed signal:
prepfigure(1);
subplot(1,2,1);
cs.xorig.colormap_freeze=false;
cs.xorig.graph_signalInverseWindowed(false);
ca=caxis;

subplot(1,2,2);
cs.x.colormap_freeze=false;
cs.x.graph_signalInverseWindowed(false);
caxis(ca);


