%% parse results
for test = 1:nrTests
    for cond = 1:nrConditions
        if isnan(results(test,cond).e_z) == 0
            % l2 errors
            e_z(test,cond) = results(test,cond).e_z;
            e_z_diag(test,cond) = results(test,cond).e_z_diag;
            e_zFast(test,cond) = results(test,cond).e_zFast;
            e_zFast_diag(test,cond) = results(test,cond).e_zFast_diag;
            e_naive(test,cond) = results(test,cond).e_naive;
            
            % l1 errors
            e_z_l1(test,cond) = results(test,cond).e_z_l1;
            e_z_diag_l1(test,cond) = results(test,cond).e_z_diag_l1;
            e_zFast_l1(test,cond) = results(test,cond).e_zFast_l1;
            e_zFast_diag_l1(test,cond) = results(test,cond).e_zFast_diag_l1;
            e_naive_l1(test,cond) = results(test,cond).e_naive_l1;
            
            % unnormalized l1 errors
            N = results(test,cond).K / results(test,cond).percSamples;
            ave_z_l1(test,cond) = results(test,cond).e_z_l1 * results(test,cond).zGT_norm1 / N;
            ave_z_diag_l1(test,cond) = results(test,cond).e_z_diag_l1 * results(test,cond).zGT_norm1 / N;
            ave_zFast_l1(test,cond) = results(test,cond).e_zFast_l1 * results(test,cond).zGT_norm1 / N;
            ave_zFast_diag_l1(test,cond) = results(test,cond).e_zFast_diag_l1 * results(test,cond).zGT_norm1 / N;
            ave_naive_l1(test,cond) = results(test,cond).e_naive_l1 * results(test,cond).zGT_norm1 / N;
                        
            % time
            cvx_time(test,cond) = results(test,cond).cvx_time;
            cvx_time_diag(test,cond) = results(test,cond).cvx_time_diag;
            timeFast(test,cond) = results(test,cond).timeFast;
            timeFast_diag(test,cond) = results(test,cond).timeFast_diag;
            naive_time(test,cond) = results(test,cond).naive_time;
        end
    end
end

isValid = [];
for i=1:size(e_z,1)
    if sum(isnan(e_z(i,:)))==0 && sum(isnan(e_naive(i,:)))==0 
       isValid = [isValid i];
    else
        warning('Nan in results')
    end
end
% l2 errors
e_z = e_z(isValid,:);
e_z_diag = e_z_diag(isValid,:);
e_zFast = e_zFast(isValid,:);
e_zFast_diag = e_zFast_diag(isValid,:);
e_naive = e_naive(isValid,:);
% l1 errors
e_z_l1 = e_z_l1(isValid,:);
e_z_diag_l1 = e_z_diag_l1(isValid,:);
e_zFast_l1 = e_zFast_l1(isValid,:);
e_zFast_diag_l1 = e_zFast_diag_l1(isValid,:);
e_naive_l1 = e_naive_l1(isValid,:);
% unnormalized l1 errors
ave_z_l1 = ave_z_l1(isValid,:);
ave_z_diag_l1 = ave_z_diag_l1(isValid,:);
ave_zFast_l1 = ave_zFast_l1(isValid,:);
ave_zFast_diag_l1 = ave_zFast_diag_l1(isValid,:);
ave_naive_l1 = ave_naive_l1(isValid,:);
% time
cvx_time = cvx_time(isValid,:);
cvx_time_diag = cvx_time_diag(isValid,:);
timeFast = timeFast(isValid,:);
timeFast_diag = timeFast_diag(isValid,:);
naive_time = naive_time(isValid,:);
      
%% plotting stuff
z_mark = ':b'; % o
z_diag_mark = ':c'; % d
zFast_mark = '-b'; % o
zFast_diag_mark = '-c' ; % d
naive_mark = '-r'; % s
dim = 24;
saveType = 'jpg'; %eps, jpg

%% PLOT ERRORS (l2 norm)
f = figure(); clf;hold on
plot(conditions,mean(e_z),z_mark,'linewidth',3)
plot(conditions,mean(e_z_diag),z_diag_mark,'linewidth',3)
plot(conditions,mean(e_zFast),zFast_mark,'linewidth',3)
plot(conditions,mean(e_zFast_diag),zFast_diag_mark,'linewidth',3)
plot(conditions,mean(e_naive),naive_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('error l2')
if isnan(e_z_diag(1))
  legend('l1-cvx','naive','l1-nesta');   
else
  legend('l1-cvx','l1-diag-cvx','l1-nesta','l1-diag-nesta', 'naive');  
end
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorl2'),saveType); 

%% PLOT ERRORS (l1 norm)
dim = 24;
f = figure(); clf;hold on
plot(conditions,mean(e_z_l1),z_mark,'linewidth',3)
plot(conditions,mean(e_z_diag_l1),z_diag_mark,'linewidth',3)
plot(conditions,mean(e_zFast_l1),zFast_mark,'linewidth',3)
plot(conditions,mean(e_zFast_diag_l1),zFast_diag_mark,'linewidth',3)
plot(conditions,mean(e_naive_l1),naive_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('error l1')
if isnan(e_z_diag_l1(1))
  legend('l1-cvx','naive','l1-nesta');   
else
  legend('l1-cvx','l1-diag-cvx','l1-nesta','l1-diag-nesta', 'naive');  
end
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorl1'),saveType); 

%% PLOT ERRORS (unnormalized l1 norm)
dim = 24;
f = figure(); clf;hold on
plot(conditions,mean(ave_z_l1),z_mark,'linewidth',3)
plot(conditions,mean(ave_z_diag_l1),z_diag_mark,'linewidth',3)
plot(conditions,mean(ave_zFast_l1),zFast_mark,'linewidth',3)
plot(conditions,mean(ave_zFast_diag_l1),zFast_diag_mark,'linewidth',3)
plot(conditions,mean(ave_naive_l1),naive_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('non-norm. error L1 [m]')
if isnan(ave_z_diag_l1(1))
  legend('l1-cvx','naive','l1-nesta');   
else
  legend('l1-cvx','l1-diag-cvx','l1-nesta','l1-diag-nesta', 'naive');  
end
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorUnnormalizedl1'),saveType); 

%% PLOT TIME
f = figure(); clf;hold on
plot(conditions,mean(cvx_time),z_mark,'linewidth',3)
plot(conditions,mean(cvx_time_diag),z_diag_mark,'linewidth',3)
plot(conditions,mean(timeFast),zFast_mark,'linewidth',3)
plot(conditions,mean(timeFast_diag),zFast_diag_mark,'linewidth',3)
plot(conditions,mean(naive_time),naive_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('time [s]') % L2
if isnan(ave_z_diag_l1(1))
  legend('l1-cvx','naive','l1-nesta');   
else
  legend('l1-cvx','l1-diag-cvx','l1-nesta','l1-diag-nesta', 'naive');  
end
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-time'),saveType); 

%% PLOT MU vs ERROR difference (unnormalized l1 norm)
dim = 24;
f = figure(); clf;
semilogx(conditions,mean(ave_zFast_l1)-mean(ave_z_l1),z_mark,'linewidth',3)
% hold on
% semilogx(conditions,mean(ave_zFast_diag_l1)-mean(ave_z_diag_l1),zFast_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('non-norm. error diff. L1 nesta-cvx [m]') % L2
legend('l1 nesta-cvx error','l1 nesta-cvx diag error');
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-errorUnnormalizedl1_vs_mu'),saveType); 


%% PLOT MU vs TIME

dim = 24;
f = figure(); clf;
semilogx(conditions,mean(timeFast),z_mark,'linewidth',3)
% hold on
% semilogx(conditions,mean(timeFast_diag),zFast_mark,'linewidth',3)
xlabel(xaxisStr); ylabel('time [s]') % L2
legend('l1-nesta','l1-diag-nesta');
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
saveas(f,horzcat(filename,'-time_vs_mu'),saveType); 

%% PLOT TIME vs ERROR difference (unnormalized l1 norm) with various MU
meantimeFast = mean(timeFast);
ave_error_diff = mean(ave_zFast_l1)-mean(ave_z_l1);

dim = 24;
f = figure(); clf;

scatter(meantimeFast,ave_error_diff,'SizeData',80,'Marker','o','MarkerFaceColor','flat')
xlabel('time [s]'); ylabel('non-norm. error L1 [m]')
grid on
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
for k = 1:nrConditions
    txt = ['\mu = ', num2str(conditions(k))];
    text(meantimeFast(k)-0.03,ave_error_diff(k)+1e-4,txt,'VerticalAlignment','bottom','fontsize',12)
end
saveas(f,horzcat(filename,'-time_vs_errorUnnormalizedl1'),saveType); 























































