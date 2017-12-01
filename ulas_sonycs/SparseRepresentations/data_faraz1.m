% compare PSNR of compressive sensing recon (Fourier->Curvelet)
% (data computed by Faraz) with hard thresholding at critical
% oversampling rate (which should be an upper limit)


% addpathContainers;
% addpathMisc;
% WavePath;
% addpathSplittingSolvers;
% setpathSonyCS( );

cvals=[4,8,16];
tab1=cell(4,8);
% data from Faraz Siddiqui, Compressive Sensing, Sony report 19-12-13
tab1(1,:)= {'barbara','cameraman','peppers','boats','harbour','goldhill','einstein','mri'};
tab1(2,:)={33.88, 31.84, 34.89, 34.70, 29.49, 32.06, 32.90, 34.10};
tab1(3,:)={27.98, 27.15, 31.34, 30.04, 25.73, 28.61, 30.51, 30.71};
tab1(4,:)={24.18, 24.88, 28.10, 27.32, 23.94, 26.94, 26.71, 27.89};
findIdx=@(str) find(~cellfun(@isempty,strfind(tab1(1,:),str)));

w=Curvelet2_clab();
w.set_wcoarsestlev(4);
w.allcurvelets=true; % best quality of true
w.fig=3;
layout=2; % 2->eachinto1figure
fig=2;
[sel,sr]=best_sterm_approx(w, layout,fig);

subplot(1,2,1);
hold on;
for j=1:length(sel)
    idx=findIdx(sel{j});
    PSNR=cell2mat(tab1(2:end,idx));
    hfig=plot(cvals,PSNR,'o');
    hold all;
end
[LEGH,OBJH,OUTH,OUTM] = legend;
htitle = get(LEGH,'Title'); 
oldtitel=get(htitle,'String');

% Add object with new handle and new legend string to legend
legh=legend([OUTH;hfig],OUTM{:},'CS Fourier \rightarrow Curvelet');
htitle = get(legh,'Title'); 
set(htitle,'String',oldtitel);



