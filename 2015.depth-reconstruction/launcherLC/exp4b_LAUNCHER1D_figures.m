%% Fangchang: generate figures for explaining wrapping and un-wrapping of signals

clear all; close all; clc

% load x and y
load ../testFM/2D_SLAM/scan261.mat
z = y;
marg = 0.5
x = x(1:5:end);
z = z(1:5:end)+.5;

%% Wrapping 
dim = 24;
f = figure(1); clf; hold on
plot(x,z,'-ok','linewidth',3,'markersize',5,'MarkerFaceColor','black');
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
ax = gca;
% margin = 0.3;
% xlim([x(1)-2*margin x(end)+margin]);
box off
grid off
% ax.XTick = [x(1) : 1 : x(end)];
xlabel('x');
yl = ylabel('z (depth)');
for i=1:length(x)
   plot([x(i) x(i)],[0 z(i)],'-b') ;
end
ylim([0 4.5])
xlim([min(x)-marg max(x)+marg])
axis equal
% title('Original Signal')

%% Un-Wrapping 
x_unwrapped = zeros(size(x));
for i = 1 : length(x)-1
    x_unwrapped(i+1) = x_unwrapped(i) + abs(x(i+1) - x(i));
end
f = figure(2); clf; hold on
plot(x_unwrapped,z,'-ok','linewidth',3,'markersize',5,'MarkerFaceColor','black');
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
ax = gca;
% margin = 0.3;
% xlim([x(1)-2*margin x(end)+margin]);
box off
grid off
% ax.XTick = [x(1) : 1 : x(end)];
xlabel('x (unwrapped)');
yl = ylabel('z (depth)');
for i=1:length(x)
   plot([x_unwrapped(i) x_unwrapped(i)],[0 z(i)],'-b') ;
end
ylim([0 4.5])
xlim([min(x_unwrapped)-marg max(x_unwrapped)+marg])
axis equal
% title('Un-wrapped Signal')

 