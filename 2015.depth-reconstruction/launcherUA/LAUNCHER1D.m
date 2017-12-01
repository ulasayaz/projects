clear all
close all
clc

nrTests = 10;
conditions = [1:15]; % nr corners
nrConditions = length(conditions);

%% noiseless
epsilon = 0;
for test = 1:nrTests
    for cond = 1:nrConditions
        
        nrCorners = conditions(cond);
        pass = false;
        while pass == false
            try
                [e_z(test,cond),e_z2(test,cond) ,e_n(test,cond),cvx_time1(test,cond),cvx_time2(test,cond),naive_time(test,cond)] =...
                    example_TV2_1D_theory_function(nrCorners,epsilon);
                pass = true;
            catch ME
                disp(ME.identifier)
                pause
            end
        end
    end
    close all
end

save('./resultsRSS/results1Dnoiseless')
dim = 24;
f = figure(1); clf;
hold on
plot(conditions,mean(e_z),'--m','linewidth',2)
plot(conditions,mean(e_z2),'-b','linewidth',2)
plot(conditions,mean(e_n),'-r','linewidth',2)
xlabel('nr. corners')
ylabel('error (%)')
set(gca,'FontSize',dim) %
ylabh=get(gca,'ylabel')
set(ylabh, 'FontSize', dim)
xlabh=get(gca,'ylabel')
set(xlabh, 'FontSize', dim)
set(f, 'paperunits', 'points' )
set(f,'papersize',[300 300]);
set(f,'PaperPositionMode','Auto')
saveas(f,'./resultsRSS/results1Dnoiseless','pdf'); 
saveas(f,'./resultsRSS/results1Dnoiseless.fig'); 