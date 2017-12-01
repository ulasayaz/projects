% Testing the reconstruction capabilities of a pair of Parseval frames
% in 1 dimension


% -------- 1D
% results: for sine ramp signal:
% --- Hadamard encoder with exponential sampling works best;
%     random sampling of standard basis gives worse results.
% --- higher daubechies are better, adding hartley as second frame
%     gives 4dB improvement w.r.t. db6 only.
% --- when combining an fft with a wavelet-decoder the rectangular
%     window ceases to be a problem (cf. fft only in TestFrameFourier1D).

if ~exist('c','var')
    c=2;
end
if ~exist('fig','var')
    fig=3;
end
if fig>0
    fig=max(3,fig);  % output consists of 3 figures
end
if fig>0
    fig=max(fig,3);
end
if ~exist('wvn','var')
    wvn='db2';
end
if ~exist('fweight1','var')
    fweight1=0.5;  % weight of frame 1
end
if ~exist('sig','var')
    sig=0; % 0.01;
end
if ~exist('use_eps','var')
    use_eps=false;
end
if ~exist('use_Hadamard','var')
    use_Hadamard=true;
end
if ~exist('eps_multiple','var')
    eps_multiple=0.25;  %0.02 for sig=0.01 and length=256
end
if ~exist('signal','var') || ~isa(signal,'TimeSignal')
    %signal=TimeSignal.make_superpos();
    signal=TimeSignal.make_sinusRamp(3,0,256);    
    signal.normalize;
    if sig>0
        signal.addNoise(sig);
    end    
end
if ~exist('params','var')
    params=struct;
end
signal.marker='';


Phi1=HartleyTrafo1(signal);
Phi2=Wavelet1D(signal); Phi2.set_basisname(wvn);


if use_Hadamard
    % exponential random sampling in Hadamard transform:
    B=WalshHadamard1(signal);  
    %nodes=B.nodesExpRnd(c);
    nodes=B.nodesDefault(c,params);    
else
    B=IdTransform(signal);
    %--- random nodes:
    %nodes=B.nodesDefault(c);
    %--- large gap in sampling:
    gap=[ceil(signal.length/3),ceil(2*signal.length/3)];    
    nodes=SampledNodes([1:gap(1),gap(2):signal.length]);
end


Phi=FrameTrafoPair(signal);
Phi.set_synweight1(fweight1);
Phi.set_transforms(Phi1,Phi2);


cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
cs.keep_reconsparse=isa(Phi, 'FrameTrafoPair');

if use_eps
    cs.paramsRecon.epsilon=eps_multiple*sig*sqrt(signal.numel);
end

cs.sim_SampledONS(Phi,B,nodes);

% draw again: cs.show_recon;

if isa(Phi, 'FrameTrafoPair');
    % show reconstruction individually for frame pair:
    cs.Phi.set_C(cs.xsparse.xn);
    cs.Phi.show_recon;
end

% compare all reconstructions:
prepfigure(fig,[],cs.figopt);
subplot(2,2,1);
cs.xorig.set_scale(PosScales.n);  % choose same scale as that of encoder
cs.xorig.graph_signal(false);
if isa(B,'IdTransform')
    hold all;     
    cs.graph_encoder(false);
    hold off;
else
    h=get(gca,'Title');
    present_titstr=get(h,'String');
    if iscellstr(present_titstr)
        present_titstr{1}=['Encoder: c=',num2str(cs.c,'%3.1f'),', ', ...
            cs.A.basisname, ' - ',present_titstr{1}];
    else
        present_titstr={['Encoder: c=',num2str(cs.c,'%3.1f'),', ', ...
            cs.A.basisname], present_titstr};
    end
    title(present_titstr,'fontsize',12);
    
end

subplot(2,2,2);
cs.graph_recon(false);

subplot(2,2,3);
% test same(!) measurement vector with only first decoder:
cs.set_Phi(Phi1);
cs.rec(); 
cs.graph_recon(false);

subplot(2,2,4);
% test same(!) measurement vector with only second decoder:
cs.set_Phi(Phi2);
cs.rec(); 
cs.graph_recon(false);






