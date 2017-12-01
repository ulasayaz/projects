close all;
clear all;
clc;
experiment = 'exp15'; %exp1, exp2, exp3, exp15
ADD_LIBRARIES

switch experiment
    case 'exp1'
        load(matFolder,'exp1/1_results1D_noiseless-N2000-nrCorn15-maxVal5-addNoise0-noiseMode_l1inf-eps0-sample_inSegments-percSamples-100-addNe1-addBo1-opt_l1inf.mat')
        xaxisStr = 'nr. corners';
        xtick = [1:2:15];
        xl = [0 16];
        yl = [-0.1 0.9];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northeast';
    case 'exp2'
        load(matFolder,'exp2/2_results1D_s2ineq_percsamples-N2000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps10-sample_uniform-percSamples100-addNe1-addBo1-opt_l1inf.mat')
        xaxisStr = 'perc. samples';
        xtick = [0:0.2:1];
        xl = [-0.05 1.05];
        yl = [-0.005 0.1];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northeast';
    case 'exp3'
        load(matFolder,'exp3/3_results1D_s2ineq_epsilon-N2000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps100-sample_uniform-percSamples5-addNe1-addBo1-opt_l1inf.mat')
        xaxisStr = 'noise level [$\varepsilon$]';
        L1str = horzcat('L1($\varepsilon$)');
        xtick = [0:0.2:1];
        xl = [-0.05 1.05];
        yl = [-0.03 0.45];
        legendStr = 'northwest';
    case 'exp15'
        load(matFolder,'exp15/15_results1D_noiseless_N_-N10000-nrCorn3-maxVal5-addNoise0-noiseMode_l1inf-eps0-sample_uniform-percSamples5-addNe1-addBo1-opt_l1inf.mat')
        xaxisStr = 'n';
        xtick = [];
        xl = [500 10500];
        yl = [0 1.1];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northwest';
    otherwise
        error('wrong experiment')
end

warning('adding time of stage 1 to stage 2')
%% parse results
for test = 1:nrTests
    for cond = 1:nrConditions 
        % unnormalized l1 errors
        N = results(test,cond).K / results(test,cond).percSamples;
        ave_zlinf_l1(test,cond) = results(test,cond).e_zlinf_l1 * results(test,cond).zGT_norm1 / N;
        ave_zlinf_s2_l1(test,cond) = results(test,cond).e_zlinf_s2_l1 * results(test,cond).zGT_norm1 / N;
        ave_zl2_l1(test,cond) = results(test,cond).e_zl2_l1 * results(test,cond).zGT_norm1 / N;
        ave_zl2_s2_l1(test,cond) = results(test,cond).e_zl2_s2_l1 * results(test,cond).zGT_norm1 / N;
        ave_n_l1(test,cond) = results(test,cond).e_n_l1 * results(test,cond).zGT_norm1 / N; 
        
        % time
        cvx_time_linf(test,cond) = results(test,cond).cvx_time_linf;
        cvx_time_linf_s2(test,cond) = results(test,cond).cvx_time_linf_s2 + cvx_time_linf(test,cond);
        cvx_time_l2(test,cond) = results(test,cond).cvx_time_l2;
        cvx_time_l2_s2(test,cond) = results(test,cond).cvx_time_l2_s2 + cvx_time_l2(test,cond);
        naive_time(test,cond) = results(test,cond).naive_time;
    end
end

zlinf_mark = '--m' %'--om'
zlinf_s2_mark = '-b' % '-b'
naive_mark = ':r' % '
dim = 26;
saveType = 'epsc'; %eps, jpg

%% unnormalized PLOT l1 error
if ~isequal(experiment,'exp15')
    f = figure(); clf; hold on
    plot(conditions,mean(ave_n_l1),naive_mark,'linewidth',3.5)
    plot(conditions,mean(ave_zlinf_l1),zlinf_mark,'linewidth',3.5)
    plot(conditions,mean(ave_zlinf_s2_l1),zlinf_s2_mark,'linewidth',3.5)
    
    if isempty(xtick)==0
        ax = gca;
        ax.XTick = xtick;
    end
    if isempty(xl)==0
        xlim(xl);
    end
    if isempty(yl)==0
        ylim(yl);
    end
    
    xlabel(xaxisStr,'Interpreter','LaTex');
    ylabel('error [m]','Interpreter','LaTex')
    leg = legend('naive',L1str,'A1','Location',legendStr);
    set(leg,'Interpreter','latex');
    set(leg,'fontsize',26);
    grid on
    set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
    set(ylabh, 'FontSize', dim); xlabh=get(gca,'xlabel');
    set(xlabh, 'FontSize', dim);
    box off
    saveas(f,horzcat(filename,'-errorUnnormalizedl1-',experiment),saveType);
end


%% PLOT TIME
if isequal(experiment,'exp15')
    f = figure(); clf;hold on
    %plot(conditions.^2,mean(naive_time),naive_mark,'linewidth',3)
    plot(conditions,mean(cvx_time_linf),zlinf_mark,'linewidth',3)
    plot(conditions,mean(cvx_time_linf_s2),zlinf_s2_mark,'linewidth',3)
    
%     if isempty(xtick)==0
%         ax = gca;
%         ax.XTick = xtick;
%     end
    if isempty(xl)==0
        xlim(xl);
    end
    if isempty(yl)==0
        ylim(yl);
    end
    
    xlabel(xaxisStr,'Interpreter','LaTex');
    ylabel('time [s]','Interpreter','LaTex') % L2
    leg = legend(L1str,'A1','Location',legendStr);
    set(leg,'Interpreter','latex');
    set(leg,'fontsize',26);
    grid on
    set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
    set(ylabh, 'FontSize', dim); xlabh=get(gca,'xlabel');
    set(xlabh, 'FontSize', dim);
    box off
    saveas(f,horzcat(filename,'-time-',experiment),saveType);
end
