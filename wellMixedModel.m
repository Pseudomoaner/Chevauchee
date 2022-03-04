clear all
% close all

v = 0.1;
lam = 0.02;
atFrac = 2/10;

tMax = 1000;

[t,y] = ode45(@(t,y)wellMixedODEs(t,y,v,lam),[0,tMax],[100*(1-atFrac),0,0,0,0,100*atFrac,0,0,0,0]);

figure
subplot(1,2,1)
ax1 = gca;
plotPopFracBars(y(1,:),ax1)
title('t = 0')

subplot(1,2,2)
ax2 = gca;
plotPopFracBars(y(end,:),ax2)
title(['t = ',num2str(tMax)])