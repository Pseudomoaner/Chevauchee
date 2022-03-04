clear all
close all

v = 0.1; %Cell velocity
lam = 0.02; %CDI firing rate
atFrac = 2/10; %Attacker fraction
rho0 = 0.0001; %Seeding density
w = 1; %Cell width

initJ = 4*sqrt(rho0)*atFrac*w;

tMax = 1000;

opts = odeset('MaxStep',1);

[t,y] = ode45(@(t,y)patchyODEs(t,y,v,lam,rho0,atFrac),[0,tMax],[100*(1-initJ),0,0,0,0,100*initJ,0,0,0,0,initJ],opts);

writePopFracBarMovie(y(:,1:10),t)

% figure
% subplot(1,2,1)
% ax1 = gca;
% plotPopFracBars(y(1,1:end-1),ax1)
% title('t = 0')
% 
% subplot(1,2,2)
% ax2 = gca;
% plotPopFracBars(y(end,1:end-1),ax2)
% title(['t = ',num2str(tMax)])