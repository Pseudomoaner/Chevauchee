clear all
close all

rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [10,58,92;31,126,193;76,178,250;156,212,252]/255;

lams = [0.02]; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
hitEfficiencies = [0.01,0.02,0.05,0.1,0.2]; %The impact of each hit on the targeted cell's ongoing growth rate
atFrac = 5.5/11; %Attacker fraction

%Dimensional parameters
dx = 10; %Granularity of coarse-grained lattice
xWidth = 300; %780 = 624 um / 0.8 um
yHeight = 300; %630 = 501 um / 0.8 um
noX = xWidth/dx;
noY = yHeight/dx;
w = 1; %Cell width
tMin = -2.5*3600;
tMax = 1.5*3600;
diffDt = 200; %Timestep between diffusion timesteps
noDiffTsteps = (tMax-tMin)/diffDt;
tList = linspace(tMin,tMax,noDiffTsteps+1);
noReps = 5;

velLists = (repmat(getExptVelocityCourses(tList/3600),1,1,noReps))/0.8;

%Rate-related parameters
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noConts = 5; %Number of contacts made by each cell

for lamInd = 1:size(lams,2)
    lam = lams(lamInd);
    for efInd = 1:size(hitEfficiencies,2)
        hitEfficiency = hitEfficiencies(efInd);
        noHitBins = ceil(1/hitEfficiency) + 1;

        subplot(size(lams,2),size(hitEfficiencies,2),(efInd-1)*size(lams,2)+lamInd)
        ax = gca;
        hold(ax,'on')

        CIstore = zeros(size(rho0s,2),noReps);
        for i = 1:4
            rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)
            
            vList = velLists(:,i);
            
            for r = 1:noReps
                [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
                pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)
                
                sensTcourse = zeros(noDiffTsteps,1);
                
                %Allow the system to equilibrate contact compartments without killing or
                %diffusion
                [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,3600],pops(:));
                pops = pops(end,:);
                pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
                
                atPopSize = sum(sum(pops(:,:,1),1),2);
                sensPopSize = sum(sum(sum(pops(:,:,2:end),1),2),3);

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
                    
                    for j = 0:noHitBins-1
                        sensTcourse(tS) = sensTcourse(tS) + squeeze(sum(sum(sum(pops(:,:,2+j:noHitBins:end),1),2),3)) * max(0,(1-j*hitEfficiency));
                    end
                end

                CIs = (atPopSize./sensTcourse)./(atPopSize/sensPopSize);
                CIstore(i,r) = CIs(end);
                
                endPops = pops;
                
                save(sprintf('C:\\Users\\Olivier\\OneDrive - Université de Lausanne\\Simulations\\Driveby\\popDat_rho_%f_ef_%f_lam_%f_r_%i.mat',rho0,hitEfficiency,lam,r),'startA','startS','CIs','endPops','sensTcourse','atPopSize','sensPopSize')
            end
        end

        plotSpread(ax,flip(CIstore',2),'DistributionMarkers','o','BinWidth',0.5,'distributionColors',flip(colours,1),'xValues',1:4)
        for i = 1:4
            plot(ax,i,mean(CIstore(4-i+1,:)),'_','MarkerSize',15,'MarkerEdgeColor',colours(4-i+1,:),'LineWidth',1.5)
        end
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        ax.YScale = 'log';
        ax.XDir = 'reverse';
        axis(ax,[0.0005,5,1,1e20])
        xticks(1:4)
        yticks([1,1e10,1e20])
        title(['\lambda = ', num2str(lam), ', e = ', num2str(hitEfficiency)])
    end
end