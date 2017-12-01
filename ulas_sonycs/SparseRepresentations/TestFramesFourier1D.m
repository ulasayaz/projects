% Testing the 1-d Fourier transform as sparsifying transform in
% compressive sensing
% RESULTS: -- fft works only if signla has been windowed by an appropriate
%             window like @chebwin.
%          -- center of signal is very well reconstructed in the examples,
%             however border of reconstruction only very badly.
%          -- without appropriate windowing fft will usually not produce a sparse
%             signal even for oscillations whose periodicity is not the
%             window size.
%         --- when combining an fft with a wavelet-decoder the rectangular
%            window ceases to be a problem (cf. TestFrameTrafoPair1D).


if ~exist('c','var')
    c=2;
end
if ~exist('fig','var')
    fig=2;
end
if ~exist('wvn','var')
    wvn='db2';
end
if~exist('do_windowing','var')
    do_windowing=true;
end
if ~exist('signalorig','var') || ~isa(signalorig,'TimeSignal')
    %signal=TimeSignal.make_superpos();
    signalorig=TimeSignal.make_sinusRamp(3,0,256);    
    signalorig.normalize;
    if exist('sig','var') && sig>0
        signal.addNoise(sig);
    end    
end
signal.marker='';

signal=signalorig.clone;

if do_windowing
    % set window to improve frequency sparsity:
    signal.set_winfun(@chebwin);
    signal=signal.apply_winfun;
end


Phi=HartleyTrafo1(signal);

% EITHER random sampling in standard basis (bad result) .....
B=IdTransform(signal);
c=5;  nodes=B.nodesDefault(c);

% ---- OR sampling in ONS, here Walsh-Hadamard
% for better result:

% B=WalshHadamard1(signal);
% nodes=B.nodesExpRnd(c);

cs=CompressiveSensing();
cs.set_original(signal);
cs.fig=fig;
cs.keep_reconsparse=isa(Phi, 'FrameTrafoPair');

cs.sim_SampledONS(Phi,B,nodes);

if  do_windowing
    % remove window weights from reconstructed signal:
    prepfigure(1, [], cs.figopt);
    cs.xorig.graph_signalInverseWindowed(false);
    ylimvals=ylim;
    hold all;
    
    cs.x.marker='';
    cs.x.graph_signalInverseWindowed(false);
    legend('original','reconstructed','location','best');
    ylim(ylimvals);
    hold off;
end
