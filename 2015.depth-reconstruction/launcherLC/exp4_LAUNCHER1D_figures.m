clear all
close all
clc

x = [-5:1:5]; x = x(:);
xq = [x(1) -1  x(end)];
yq = [1 3  1.5]-0.5;
% generate linear interpolation
z = interp1(xq,yq,x,'linear','extrap')'; z = z(:);

dim = 24;
f = figure(); clf; hold on
plot(x,z,'-ok','linewidth',3,'markersize',7,'MarkerFaceColor','black');
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
ax = gca;
ylim([0 3])
margin = 0.3;
xlim([x(1)-2*margin x(end)+margin]);
box off
grid off
ax.XTick = [x(1) : 1 : x(end)];
xlabel('x');
yl = ylabel('z (depth)');
for i=1:length(x)
   plot([x(i) x(i)],[0 z(i)],'-b') ;
end


[xpol,r] = cart2pol(x,z);
xpol = -xpol;
dim = 24;
f = figure(); clf; hold on
plot(xpol,z,'-or','linewidth',3,'markersize',7,'MarkerFaceColor','black');
plot(xpol,r,'--om','linewidth',3,'markersize',7,'MarkerFaceColor','black');
set(gca,'FontSize',dim); ylabh=get(gca,'ylabel');
set(ylabh, 'FontSize', dim); xlabh=get(gca,'ylabel');
set(xlabh, 'FontSize', dim);
ax = gca;
ax.XTick = [-2:0.1:2];
ylim([0 5.5])
margin = 0.3;
xlim([min([xpol(1) xpol(end)])-margin  max([xpol(1) xpol(end)])+margin]);
box off
grid off
ax.XTick = [-pi -pi/2 0];
ax.XTickLabel = {'\pi','\pi/2','0'};
xlabel('\theta');
%yl = ylabel('z (depth)');
leg = legend('z','r')
set(leg,'fontsize',26);
for i=1:length(x)
   plot([xpol(i) xpol(i)],[0 r(i)],'-b');
end
 