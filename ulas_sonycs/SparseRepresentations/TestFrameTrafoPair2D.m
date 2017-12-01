% Testing the reconstruction capabilities of a pair of Parseval frames
% in 2 dimension2.
% Example:
% ----------
% 1) define signal:   signal=Signal2D.make_fromImage('cameraman.bmp');
% 2) define decoder 1:  Phi1=Curvelet2_clab(signal);
%               OR: Phi1=Fourier2(signal);
% 3) take default decoder 2 or define one yourself.
% 4) fix analysis weight for Phi1: aweight1=0.5;  
% 5) take default encoder or define one yourself: B=Wavelet2D_mlab(signal);
%              (needed for wavelet encoder only): B.set_deepestlev(1);
%                                                 params.w1=0.999;
%

if ~exist('c','var')
    c=8;  % measurement compression (NOT sparisfying compression)
end
if ~exist('RL','var')
    RL=50; 
end

if ~exist('fig','var')
    fig=3;
end
if fig>0
    fig=max(fig,3);  % output needs at least 3 images
end

if ~exist('wvn','var')
    wvn1='db5';
end
if ~exist('wvn','var')
    wvn2='db1';
end
if ~exist('params','var')
    params=struct;
end

if ~exist('signal','var') || ~isa(signal,'Signal2D')   
    signal=Signal2D.make_fromImage('cameraman.bmp');
end



if ~exist('B','var')
    %B=DCosineTrafo2(signal);
    B=WalshHadamard2(signal);  
end
if ~exist('anweight1','var')
    % determines influence of each transform
    anweight1=0.5;
end
if ~exist('Phi1','var')
    Phi1=Wavelet2D_mlab(signal); 
    Phi1.set_basisname(wvn1);
end
if ~exist('Phi2','var')
    Phi2=Wavelet2D_mlab(signal);
    Phi2.set_basisname(wvn2);
    %Phi2=Curvelet2_clab(signal);
end

nodes=B.nodesDefault(c,params);

Phi=FrameTrafoPair(signal);
Phi.set_transforms(Phi1,Phi2);
Phi.set_anweight1(anweight1);  % determines influence of each transform

cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
% set parameter to store sparsified signal in cs.xsparse:
cs.keep_reconsparse=isa(Phi, 'FrameTrafoPair');

cs.sim_SampledONS(Phi,B,nodes);

% draw again: cs.show_recon;

if isa(Phi, 'FrameTrafoPair');
    % show reconstruction individually for frame pair:
    cs.Phi.set_C(cs.xsparse.xn);
    cs.Phi.show_recon;
end

% test same(!) measurement vector with only second decoder:
cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
cs.sim_SampledONS(Phi2,B,nodes);     

% test same(!) measurement vector with only first decoder:
cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
cs.sim_SampledONS(Phi1,B,nodes);     




