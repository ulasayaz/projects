close all
clear all
clc
experiment = 'exp5'; % exp5d, exp6b, exp10b, exp11b, exp16, exp17, exp19b, exp19c, exp20b
ADD_LIBRARIES

switch experiment
   case 'exp5'
        load(matFolder,'exp5/5_results2D_newNoiseModel_percsamples-N10000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps10-sample_uniform-percSamples100-addNe0-addBo0-opt_l1inf-diag1-DIAG05.mat')
        xaxisStr = 'perc. samples';
        xtick = [0:0.2:1];
        ytick = [0:0.02:0.12];
        xl = [-0.05 1.05];
        yl = [0 0.11];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northeast';
    case 'exp5d'
        load(matFolder,'exp5d/5d_results2D_newNoiseModel_percsamples-N10000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps10-sample_uniform-percSamples100-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'perc. samples';
        xtick = [0:0.2:1];
        xl = [-0.05 1.05];
        yl = [0 0.08];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northeast';
    case 'exp6b'
        load(matFolder,'exp6b/6b_results2D_newNoiseModel_epsilon-N10000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps100-sample_uniform-percSamples5-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'noise level [$\varepsilon$]';
        xtick = [0:0.2:1];
        xl = [-0.05 1.05];
        yl = [-0.03 0.45];
        L1str = horzcat('L1($\varepsilon$)');
        L1diagstr = horzcat('L1diag($\varepsilon$)');
        legendStr = 'northwest';
    case 'exp10b'
        load(matFolder,'exp10b/10b_results2D_gazebo_percSamples-N10000-nrCorn-1-maxVal-1-addNoise0-noiseMode_l1inf-eps0-sample_uniform-percSamples100-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'perc. samples';
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        xtick = [0:0.1:0.4];
        xl = [-0.02 0.42];
        yl = [-0.01 0.22];
        legendStr = 'northeast';
        conditions = [0.02:0.02:0.4];
        nrConditions = numel(conditions);
    case 'exp11b'
        load(matFolder,'exp11b/11b_results2D_zed_percSamples-N10000-nrCorn-1-maxVal-1-addNoise0-noiseMode_l1inf-eps0-sample_uniform-percSamples100-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'perc. samples';
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        xtick = [0:0.1:0.4];
        xl = [-0 0.42];
        yl = [-0.01 0.16];
        legendStr = 'northeast';
        conditions = [0.02:0.02:0.4];
        nrConditions = numel(conditions);
    case 'exp16'
        load(matFolder,'exp16-10tests/16_results2D_noiseless_N-N90000-nrCorn3-maxVal5-addNoise0-noiseMode_l1inf-eps0-sample_uniform-percSamples5-addNe0-addBo0-opt_l1inf-diag1-DIAG05.mat')
        xaxisStr = 'n';
        %xtick = [50:50:90000];
        xl = [40 93000];
        yl = [-10 120];
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        legendStr = 'northwest';
    case 'exp17'
        load(matFolder,'exp17/17_results2D_temporal_horizon-N10000-nrCorn3-maxVal5-addNoise1-noiseMode_l1inf-eps10-sample_uniform-percSamples5-addNe0-addBo0-opt_l1inf-diag1-DIAG05.mat')
        xaxisStr = 'horizon';
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        xtick = [0:5:20];
        xl = [-0 21];
        yl = [0.04 0.13];
        legendStr = 'northwest';
    case 'exp19b'
        load(matFolder,'exp19b/19b_results2D_gazeboOrtho_percSamples-N10000-nrCorn-1-maxVal-1-addNoise0-noiseMode_l1inf-eps0-sample_uniform-percSamples100-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'perc. samples';
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        xtick = [0:0.1:0.3];
        xl = [-0 0.32];
        yl = [-0.01 0.23];
        legendStr = 'northeast';
        conditions = [0.02:0.02:0.3];
        nrConditions = numel(conditions);
    case 'exp19c'
        load(matFolder,'exp19c/19c_results2D_gazeboOrtho_percSamples_noisy-N10000-nrCorn-1-maxVal-1-addNoise1-noiseMode_l1inf-eps10-sample_uniform-percSamples100-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'perc. samples';
        L1str = horzcat('L1 ($\varepsilon$=',num2str(settings.epsilon),')');
        L1diagstr = horzcat('L1diag ($\varepsilon$=',num2str(settings.epsilon),')');
        xtick = [-0:0.2:1];
        xl = [-0.05 1.05];
        yl = [-0.01 0.23];
        legendStr = 'northeast';
    case 'exp20b'
        load(matFolder,'exp20b/20b_results2D_gazeboOrtho_epsilon-N10000-nrCorn-1-maxVal-1-addNoise1-noiseMode_l1inf-eps94-sample_uniform-percSamples5-addNe1-addBo0-opt_l1inf-diag1-DIAG05-noEff.mat')
        xaxisStr = 'noise level [$\varepsilon$]';
        L1str = horzcat('L1($\varepsilon$)');
        L1diagstr = horzcat('L1diag($\varepsilon$)');
        xtick = [0:0.2:1];
        xl = [-0.05 1.05];
        yl = [-0 0.5];
        legendStr = 'northwest';
    otherwise
        error('wrong experiment')
end

%% parse results
for test = 1:nrTests
    for cond = 1:nrConditions
        if isnan(results(test,cond).e_z) == 0
            % l2 errors
            e_z(test,cond) = results(test,cond).e_z;
            e_z_diag(test,cond) = results(test,cond).e_z_diag;
            e_n(test,cond) = results(test,cond).e_n;
            
            % unnormalized l1 errors
            N = results(test,cond).K / results(test,cond).percSamples;
            ave_z_l1(test,cond) = results(test,cond).e_z_l1 * results(test,cond).zGT_norm1 / N;
            ave_z_diag_l1(test,cond) = results(test,cond).e_z_diag_l1 * results(test,cond).zGT_norm1 / N;
            ave_n_l1(test,cond) = results(test,cond).e_n_l1 * results(test,cond).zGT_norm1 / N;
            
            % time
            cvx_time1(test,cond) = results(test,cond).cvx_time1;
            cvx_time_diag(test,cond) = results(test,cond).cvx_time_diag;
            naive_time(test,cond) = results(test,cond).naive_time;
        end
    end
end

isValid = [];
for i=1:size(ave_z_diag_l1,1)
    if sum(isnan(ave_z_diag_l1(i,:)))==0 && sum(isnan(ave_n_l1(i,:)))==0 && sum(isnan(ave_z_l1(i,:)))==0
        isValid = [isValid i];
    else
        warning('Nan in results')
    end
end

ave_z_l1 = ave_z_l1(isValid,:);
ave_z_diag_l1 = ave_z_diag_l1(isValid,:);
ave_n_l1 = ave_n_l1(isValid,:);
% time
cvx_time1 = cvx_time1(isValid,:);
cvx_time_diag = cvx_time_diag(isValid,:);
naive_time = naive_time(isValid,:);

%% plotting stuff
z_mark = '--m'
z_diag_mark = '-b'
naive_mark = ':r'
dim = 26;
saveType = 'epsc'; %eps, jpg

%% PLOT ERRORS (unnormalized l1 norm)
if ~isequal(experiment,'exp16')
    f = figure(); clf; hold on
    plot(conditions,mean(ave_n_l1),naive_mark,'linewidth',3.5)
    plot(conditions,mean(ave_z_l1),z_mark,'linewidth',3.5)
    plot(conditions,mean(ave_z_diag_l1),z_diag_mark,'linewidth',3.5)
    
    if isempty(xtick)==0
        ax = gca;
        ax.XTick = xtick;
    end
    if isempty(ytick)==0
        ax = gca;
        ax.YTick = ytick;
    end
    if isempty(xl)==0
        xlim(xl);
    end
    if isempty(yl)==0
        ylim(yl);
    end
    
    xlabel(xaxisStr,'Interpreter','LaTex');
    ylabel('error [m]','Interpreter','LaTex')
    leg = legend('naive',L1str,L1diagstr,'Location',legendStr);
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
if isequal(experiment,'exp16')
    f = figure(); clf;hold on
    %plot(conditions.^2,mean(naive_time),naive_mark,'linewidth',3)
    plot(conditions.^2,mean(cvx_time1),z_mark,'linewidth',3)
    plot(conditions.^2,mean(cvx_time_diag),z_diag_mark,'linewidth',3)
    
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
    leg = legend(L1str,L1diagstr,'Location',legendStr);
    set(leg,'Interpreter','latex');
    set(leg,'fontsize',26);
    grid on
    set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
    set(ylabh, 'FontSize', dim); xlabh=get(gca,'xlabel');
    set(xlabh, 'FontSize', dim);
    box off
    saveas(f,horzcat(filename,'-time-',experiment),saveType);
end

