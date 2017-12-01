% load('C:\GT\Projekte\2011-10_Zeiss-SMT-MIK&MO\Matlab\Pruefturm\Ergebnisse\delta2,1.mat');
load('C:\GT\Projekte\2011-10_Zeiss-SMT-MIK&MO\Matlab\Pruefturm\Ergebnisse\calcDPL5_withdpl7.mat');
h=resultDPL.phi(:,:,3)-resultDPL.phi(:,:,2);

signal=Signal2D(h);
signal.outliers2NaN(1e-3);
signal.signalname='Schliere1_{3,2}';
signal.graph_signal;
signal.hide_outliers=true;
% save('schlieren1_23.mat','signal');

nodes=signal.signal2nodes;
measvec=signal.xn(nodes);
signalsize=signal.size;

B=IdTransform();
Phi1=HartleyTrafo2();
Phi2=Wavelet2D_mlab();
wvn='db3';
Phi2.set_basisname(wvn); % 'db1' ist glatt bis auf wenige Sprünge, 'db2'

% Phi=FrameTrafoPair();
% Phi.set_transforms(Phi1,Phi2);

Phi=Phi2;

cs=CompressiveSensing();
so=Signal2D();
so.colormap_active='default';
so.signalname=signal.signalname;
so.hide_outliers=true;
cs.set_original(so);

sig=0.01;
eps= 0.25*sig*sqrt(length(nodes));
cs.paramsRecon.epsilon=eps;
cs.fig=2;

cs.compute_SampledONS(Phi,B,nodes,measvec,signalsize);
% cs.x.hide_outliers=true;
% cs.show_recon();

% unfreezeColors
% caxis([0.35,0.62]);

% cs.set_APhiNodes(Phi,B,nodes,signalsize);
% cs.set_y(measvec(:));
% cs.graph_encoder();

