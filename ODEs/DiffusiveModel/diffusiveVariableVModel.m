clear all
close all

vmaxes = [0.005,0.01,0.035,0.032]/0.8;
confTs = [4.5,9.5,11.5,14]*3600;
rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [0,0,0;0.5,0.5,0.5;0.75,0.75,0.75;0.875,0.875,0.875];

lam = 0.005; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
atFrac = 5.5/11; %Attacker fraction
vrate = 3600; %Width of Gaussian velocity profile (s^-1)

%Dimensional parameters
dx = 10; %Granularity of coarse-grained lattice
xWidth = 780; % =624 um / 0.8 um
yHeight = 630; % = 501 um / 0.8 um
noX = xWidth/dx;
noY = yHeight/dx;
w = 1; %Cell width
tMax = 24*3600;
diffDt = 200; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt;
tList = linspace(0,tMax,noDiffTsteps+1);

%Rate-related parameters
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

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

    [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
    pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)
       
    %Allow the system to equilibrate contact compartments without killing or
    %diffusion
    [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,3600],pops(:));
    pops = pops(end,:);
    pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
    popsTcourse = pops;

    for t = 1:noDiffTsteps
        currV = vmax*exp(-((t*diffDt-confTime)/vrate)^2);
        currD = currV*alphaD;

        %Run diffusion of each of the populations separately
        for j = 1:size(pops,3)
            pops(:,:,j) = diffTimestepCN(pops(:,:,j),diffDt,dx,currD,true);
        end

        %Inner loop - microscopic mixing (contact swapping)
        [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,currV,lam,noX,noY,noConts,noHitBins),[0,diffDt],pops(:));
        pops = pops(end,:);
        pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
        popsTcourse = cat(4,popsTcourse,pops);
    end
    
    atPopSize = sum(sum(popsTcourse(:,:,1,1),1),2);
    sensPopSize = sum(sum(sum(popsTcourse(:,:,2:end,1),1),2),3);

    plot(ax1,tList/3600,vmax*exp(-((tList-confTime)/vrate).^2),'LineWidth',1.5,'Color',colours(i,:))

    J = squeeze(popsTcourse(:,:,1,:));
    plot(ax2,tList/3600,squeeze(mean(mean(J.*(1-J),1),2))./(sensPopSize/(atPopSize+sensPopSize)),'LineWidth',1.5,'Color',colours(i,:))
    
    unhitTcourse = squeeze(sum(sum(sum(popsTcourse(:,:,2:noHitBins:end,:),1),2),3));
    plot(ax3,tList/3600,(atPopSize./unhitTcourse)./(atPopSize/sensPopSize),'LineWidth',1.5,'Color',colours(i,:))

    save(sprintf('C:\\Users\\olijm\\Desktop\\SeanAna\\popTcourse_%i.mat',i),'popsTcourse','vmax','rho0','confTime')
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

axis(ax1,[0,24,0,0.05])
axis(ax2,[0,24,0,0.55])
axis(ax3,[0,24,0.5,1e12])

ax3.YScale = 'log';

legend(ax1,'OD 1','OD 0.1','OD 0.01','OD 0.001')