clear all
close all

vmaxes = [0.005,0.016,0.029,0.032]/0.8;
confTs = [7200,7200,7200,7200];%[4.5,7.5,10.5,14]*3600;
rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [10,58,92;31,126,193;76,178,250;156,212,252]/255;

lams = [0.0003,0.001,0.003,0.01,0.03]; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
hitEfficiencies = [0.02,0.05,0.1,0.2,0.5]; %The impact of each hit on the targeted cell's ongoing growth rate
atFrac = 5.5/11; %Attacker fraction
vrate = 3000; %Width of Gaussian velocity profile (s^-1)

%Dimensional parameters
dx = 10; %Granularity of coarse-grained lattice
xWidth = 340; %780 = 624 um / 0.8 um
yHeight = 340; %630 = 501 um / 0.8 um
noX = xWidth/dx;
noY = yHeight/dx;
w = 1; %Cell width
tMin = -2.5*3600;
tMax = 2.5*3600;
diffDt = 200; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt;
tList = linspace(0,tMax,noDiffTsteps+1);
noReps = 1;

%Option 1: Use raw velocity timecourses
% velLists = (getExptVelocityCourses(tList/3600)-0.0018)/0.8;

%Option 2: Use averaged velocity timecourses
velLists = (repmat(mean(getExptVelocityCourses(tList/3600),3),1,1,noReps)-0.0018)/0.8;

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
        for i = 1:4
            rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)
            
            vList = velLists(:,i);

            [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
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

            sensTcourse = zeros(size(popsTcourse,4),1);
            for j = 0:noHitBins-1
                sensTcourse = sensTcourse + squeeze(sum(sum(sum(popsTcourse(:,:,2+j:noHitBins:end,:),1),2),3)) * max(0,(1-j*hitEfficiency));
            end

            atPopSize = sum(sum(popsTcourse(:,:,1,1),1),2);
            sensPopSize = sum(sum(sum(popsTcourse(:,:,2:end,1),1),2),3);
            plot(ax,rho0/0.064,(atPopSize./sensTcourse(end))./(atPopSize/sensPopSize),'.','Color',colours(i,:),'MarkerSize',15)

            %     save(sprintf('C:\\Users\\olijm\\Desktop\\SeanAna\\popTcourse_%i.mat',i),'popsTcourse','vmax','rho0','confTime')
        end
        ax.LineWidth = 1.5;
        ax.Box = 'on';
        ax.YScale = 'log';
        ax.XScale = 'log';
        axis(ax,[0.0005,5,1,1e8])
        title(['\lambda = ', num2str(lam), ', e = ', num2str(hitEfficiency)])
    end
end