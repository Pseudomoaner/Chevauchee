clear all
% close all

v = 0.1;
lam = 0.01;
atFrac = 1/10;

noHitBins = 6;
noConts = 5;

tMax = 500;

[t,y] = ode45(@(t,y)wellMixedODEs(t,y,v,lam,noConts,noHitBins),[0,tMax],[atFrac,getStartingPopulations(atFrac,1,noHitBins,noConts)']);

writePopFracBarMovie(y,t,noConts,noHitBins)

% figure
% subplot(1,2,1)
% ax1 = gca;
% plotPopFracBars(y(1,:),ax1)
% title('t = 0')
% 
% subplot(1,2,2)
% ax2 = gca;
% plotPopFracBars(y(end,:),ax2)
% title(['t = ',num2str(tMax)])