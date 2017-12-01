warning('adding time of stage 1 to stage 2')
%% parse results
for test = 1:nrTests
    for cond = 1:nrConditions
        % l2 errors
        e_zlinf(test,cond) = results(test,cond).e_zlinf;
        e_zlinf_s2(test,cond) = results(test,cond).e_zlinf_s2;
        e_zl2(test,cond) = results(test,cond).e_zl2;
        e_zl2_s2(test,cond) = results(test,cond).e_zl2_s2;
        e_n(test,cond) = results(test,cond).e_n;
        
        % l1 errors
        e_zlinf_l1(test,cond) = results(test,cond).e_zlinf_l1;
        e_zlinf_s2_l1(test,cond) = results(test,cond).e_zlinf_s2_l1;
        e_zl2_l1(test,cond) = results(test,cond).e_zl2_l1;
        e_zl2_s2_l1(test,cond) = results(test,cond).e_zl2_s2_l1;
        e_n_l1(test,cond) = results(test,cond).e_n_l1;
 
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


zlinf_mark = '--dc'
zlinf_s2_mark = '-db'
zl2_mark = '--ok'
zl2_s2_mark = '-om'
naive_mark = '-sr'
dim = 24;
saveType = 'jpg'; %eps, jpg
%% PLOT l2 error
f = figure(); clf; hold on
plot(conditions,mean(e_zlinf),zlinf_mark,'linewidth',3)
plot(conditions,mean(e_zlinf_s2),zlinf_s2_mark,'linewidth',3)
plot(conditions,mean(e_zl2),zl2_mark,'linewidth',2)
plot(conditions,mean(e_zl2_s2),zl2_s2_mark,'linewidth',2)
plot(conditions,mean(e_n),naive_mark,'linewidth',2)
xlabel(xaxisStr); ylabel('error (L2)')
legend('linf','linf-s2','l2','l2-s2','naive'); grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorl2'),saveType); 

%% PLOT l1 error
f = figure(); clf; hold on
plot(conditions,mean(e_zlinf_l1),zlinf_mark,'linewidth',3)
plot(conditions,mean(e_zlinf_s2_l1),zlinf_s2_mark,'linewidth',3)
plot(conditions,mean(e_zl2_l1),zl2_mark,'linewidth',2)
plot(conditions,mean(e_zl2_s2_l1),zl2_s2_mark,'linewidth',2)
plot(conditions,mean(e_n_l1),naive_mark,'linewidth',2)
xlabel(xaxisStr); ylabel('error (L1)')
legend('linf','linf-s2','l2','l2-s2','naive'); grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorl1'),saveType); 

%% unnormalized PLOT l1 error
f = figure(); clf; hold on
plot(conditions,mean(ave_zlinf_l1),zlinf_mark,'linewidth',3)
plot(conditions,mean(ave_zlinf_s2_l1),zlinf_s2_mark,'linewidth',3)
plot(conditions,mean(ave_zl2_l1),zl2_mark,'linewidth',2)
plot(conditions,mean(ave_zl2_s2_l1),zl2_s2_mark,'linewidth',2)
plot(conditions,mean(ave_n_l1),naive_mark,'linewidth',2)
xlabel(xaxisStr); ylabel('non-norm. error L1 [m]')
legend('linf','linf-s2','l2','l2-s2','naive'); grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorUnnormalizedl1'),saveType); 

%% PLOT TIME
f = figure(); clf;hold on
plot(conditions,mean(cvx_time_linf),zlinf_mark,'linewidth',3)
plot(conditions,mean(cvx_time_linf_s2),zlinf_s2_mark,'linewidth',3)
plot(conditions,mean(cvx_time_l2),zl2_mark,'linewidth',3)
plot(conditions,mean(cvx_time_l2_s2),zl2_s2_mark,'linewidth',3)
plot(conditions,mean(naive_time),naive_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('time [s]') % L2
legend('linf','linf-s2','l2','l2-s2','naive');
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-time'),saveType); 

% f = figure(); clf;hold on
% plot(conditions, ( mean(e_n) - mean(e_zlinf) ) ./ mean(e_n),'--m','linewidth',2)
% plot(conditions, ( mean(e_n) - mean(e_zlinf_s2)) ./ mean(e_n),zlinf_s2_mark,'linewidth',2)
% xlabel(xaxisStr); ylabel('improv. wrt naive (L2)')
% legend('L1','alg1'); grid on
% set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
% set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
% set(xlabh, 'FontSize', dim);

% f = figure(); 
% semilogy(conditions,mean(e_zlinf),'--m','linewidth',2); hold on
% semilogy(conditions,mean(e_zlinf_s2),zlinf_s2_mark,'linewidth',2)
% semilogy(conditions,mean(e_n),naive_mark,'linewidth',2)
% xlabel(xaxisStr); ylabel('error (L2)')
% legend('L1','alg1','naive'); grid on
% set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
% set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
% set(xlabh, 'FontSize', dim);