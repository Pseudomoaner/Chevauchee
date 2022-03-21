clear all
close all

vmaxes = [0.003,0.006,0.01,0.01]/0.8;
confTs = [4.5,9.5,11.5,14]*3600;
rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [1,0.8,0;1,0.2,0;0,0.5,1;0,0,1];

lam = 0.01; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
atFrac = 5.5/11; %Attacker fraction
vrate = 3600; %Width of Gaussian velocity profile (s^-1)
totPopSize = 100;
atPopSize = totPopSize*atFrac;
sensPopSize = totPopSize-atPopSize;

figure(1);
subplot(3,1,1)
ax1=gca;
hold on
subplot(3,1,2)
ax2=gca;
hold on
subplot(3,1,3)
ax3=gca;
hold on

for i = 1:4
    rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)
    vmax = vmaxes(i); %Maximum velocity (cewllWidths s^-1)
    confTime = confTs(i); %Time to confluency

    w = 1; %Cell width
    noHitBins = 5;

    initJ = 4*sqrt(rho0)*atFrac*w;

    tMax = 24*3600;

    opts = odeset('MaxStep',1);

    [t,y] = ode45(@(t,y)patchyODEsVariableV(t,y,vmax,vrate,confTime,lam,rho0,atFrac),[0,tMax],[getStartingPopulations(initJ,sensPopSize,noHitBins)',initJ],opts);

    plot(ax1,t/3600,vmax*exp(-((t-confTime)/vrate).^2),'LineWidth',1.5,'Color',colours(i,:))
    plot(ax2,t/3600,y(:,end),'LineWidth',1.5,'Color',colours(i,:))
    plot(ax3,t/3600,(atPopSize./sum(y(:,1:5:21),2))./(atPopSize/sensPopSize),'LineWidth',1.5,'Color',colours(i,:))
end

xlabel(ax1,'Time (hr)')
ylabel(ax1,'Velocity (cell widths s^{-1})')
xlabel(ax2,'Time (hr)')
ylabel(ax2,'J')
xlabel(ax3,'Time (hr)')
ylabel(ax3,{'Competitve advantage of CDI+'})

ax1.Box = 'on';
ax1.LineWidth = 1.5;
ax2.Box = 'on';
ax2.LineWidth = 1.5;
ax3.Box = 'on';
ax3.LineWidth = 1.5;

axis(ax1,[0,24,0,0.014])
axis(ax2,[0,24,0,0.11])
axis(ax3,[0,24,0.5,1e5])

ax3.YScale = 'log';

legend(ax1,'OD 1','OD 0.1','OD 0.01','OD 0.001')