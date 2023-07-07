clear all
close all

rng(0) %Allows reproducible generation of the initial starting configurations

rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [10,58,92;31,126,193;76,178,250;156,212,252]/255;
coloursMixes = [198, 159, 137;206, 45, 45;119, 225, 119]/255;

lam = 0.005; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
hitEfficiency = 0.02; %The impact of each hit on the targeted cell's ongoing growth rate
atFrac = 5.5/11; %Attacker fraction

%Dimensional parameters
dx = 10; %Granularity of coarse-grained lattice
xWidth = 350; %780 = 624 um / 0.8 um
yHeight = 350; %630 = 501 um / 0.8 um
noX = xWidth/dx;
noY = yHeight/dx;
w = 1; %Cell width
tMin = -2.5*3600;
tMax = 1.5*3600;
diffDt = 200; %Timestep between diffusion timesteps
noDiffTsteps = (tMax-tMin)/diffDt;
tList = linspace(tMin,tMax,noDiffTsteps+1);
noReps = 5; %Must be set to three if the raw velocity timecourses are used (see next line)

%Option 1: Use raw velocity timecourses
% velLists = (getExptVelocityCourses(tList/3600))/0.8;

%Option 2: Use averaged velocity timecourses
velLists = (repmat(mean(getExptVelocityCourses(tList/3600),3),1,1,noReps))/0.8;

%Rate-related parameters
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noHitBins = ceil(1/hitEfficiency) + 1; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

figure(1);
subplot(2,1,1)
ax1=gca;
hold on
subplot(2,1,2)
ax2=gca;
hold on

figure(2);
ax21 = gca;
hold on

figure(3);
ax31 = gca;

baseStore = zeros(size(rho0s,2),noReps);
contExStore = zeros(size(rho0s,2),noReps);
homoStore = zeros(size(rho0s,2),noReps);

CIstore = zeros(size(rho0s,2),noReps);

Jstore = zeros(noDiffTsteps+1,size(rho0s,2),noReps);

for i = 1:4
    for r = 1:noReps
        rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)

        vList = velLists(:,i,r);

        [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);

        [baseStore(i,r),contExStore(i,r),homoStore(i,r)] = calcFullMixingContributions(dx,startA,startS,atFrac,noConts,noHitBins,noX,noY,lam,vList,diffDt,tMax-tMin,hitEfficiency);

        pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

        %Allow the system to equilibrate contact compartments without killing or
        %diffusion
        [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,3600],pops(:));
        pops = pops(end,:);
        pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
        popsTcourse = pops;

        for tS = 1:noDiffTsteps
            currV = vList(tS);
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

        yyaxis(ax1,'left')
        plot(ax1,tList/3600,vList*alphaD,'-','LineWidth',1.5,'Color',colours(i,:))
        yyaxis(ax1,'right')
        plot(ax1,tList/3600,vList*0.1033,'-','LineWidth',1.5,'Color',colours(i,:))

        J = squeeze(popsTcourse(:,:,1,:));
        plot(ax2,tList/3600,1-var(reshape(J,noX*noY,size(J,3)))*4,'LineWidth',0.5,'Color',1-((1-colours(i,:))/2))
        Jstore(:,i,r) = 1-var(reshape(J,noX*noY,size(J,3)))*4;

        sensTcourse = zeros(size(popsTcourse,4),1);
        for j = 0:noHitBins-1
            sensTcourse = sensTcourse + squeeze(sum(sum(sum(popsTcourse(:,:,2+j:noHitBins:end,:),1),2),3)) * max(0,(1-j*hitEfficiency));
        end
        CIs = (atPopSize./sensTcourse)./(atPopSize/sensPopSize);

        CIstore(i,r) = CIs(end);
    end
end

%Figure 1 stuff
xlabel(ax1,'Time (hr)')
yyaxis(ax1,'left')
ylabel(ax1,'D')
yyaxis(ax1,'right')
ylabel(ax1,'r_e')
xlabel(ax2,'Time (hr)')
ylabel(ax2,'Strain intermixing')

ax1.Box = 'on';
ax1.LineWidth = 1.5;
ax2.Box = 'on';
ax2.LineWidth = 1.5;

yyaxis(ax1,'left')
axis(ax1,[tMin/3600,tMax/3600,0,alphaD*0.06])
yyaxis(ax1,'right')
axis(ax1,[tMin/3600,tMax/3600,0,0.1033*0.06])

axis(ax2,[tMin/3600,tMax/3600,0,1])

for i = 1:4
    plot(ax2,tList/3600,mean(Jstore(:,i,:),3),'Color',colours(i,:),'LineWidth',2.5)
end

legend(ax1,'OD 1','OD 0.1','OD 0.01','OD 0.001')

%Figure 2 stuff
meanBase = mean(baseStore,2);
meanContEx = mean(contExStore,2);
meanHomo = mean(homoStore,2);

br = bar(ax21,[meanBase,meanContEx,meanHomo],'stacked','LineWidth',1);
for k = 1:3
    br(k).FaceColor = coloursMixes(k,:);
end

plot(ax21,[0,numel(rho0s)+1],1-[atFrac,atFrac],'k--','LineWidth',1.5)
for i = 1:4
    plotSpread(ax21,baseStore','DistributionMarkers','o','distributionColors','k','BinWidth',0.5)
    plotSpread(ax21,baseStore' + contExStore','DistributionMarkers','o','distributionColors','k','BinWidth',0.5)
end

ax21.XTickLabel = cellstr(num2str(rho0s(:)/0.064));
axis(ax21,[0.4,4.6,0,1])
xlabel(ax21,'Seeding density')
ylabel(ax21,'Contribution to killing')
ax21.Box = 'on';
ax21.LineWidth = 1.5;
legend(ax21,'Base','Contact exchange','Homogenization')
xlabel(ax21,'Inoculation density')
ylabel(ax21,'Fractional contribution to killing')

%Figure 3 stuff
plotSpread(ax31,CIstore','DistributionMarkers','o','BinWidth',0.5,'distributionColors',colours)
for i = 1:4
    plot(ax31,i,mean(CIstore(i,:)),'_','MarkerSize',15,'MarkerEdgeColor',colours(i,:),'LineWidth',1.5)
end
ax31.YScale = 'log';
ax31.Box = 'on';
ax31.LineWidth = 1.5;
axis(ax31,[0.5,4.5,1,1e7])
ax31.XTickLabel = string(rho0s/0.064);
xlabel(ax31,'Inoculation density')
ylabel(ax31,'Predicted comptitive index')

re = velLists(:,:,1)*0.1033;
D = velLists(:,:,1)*alphaD;

%Save stuff
save('C:\\Users\\Olivier\\OneDrive - Université de Lausanne\\Simulations\\Driveby\\Fig4Sims.mat','Jstore','CIstore','baseStore','contExStore','homoStore','re','D')