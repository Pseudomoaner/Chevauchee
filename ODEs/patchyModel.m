clear all
close all

v = 0.1; %Cell velocity
lam = 0.01; %CDI firing rate
atFrac = 3/10; %Attacker fraction
rho0 = 0.001; %Seeding density
w = 1; %Cell width
noHitBins = 6;

initJ = 4*sqrt(rho0)*atFrac*w;

tMax = 500;

opts = odeset('MaxStep',1);

[t,y] = ode45(@(t,y)patchyODEs(t,y,v,lam,rho0,atFrac),[0,tMax],[getStartingPopulations(initJ,100,noHitBins)',initJ],opts);

writePopFracBarMovie(y(:,1:end-1),t)

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