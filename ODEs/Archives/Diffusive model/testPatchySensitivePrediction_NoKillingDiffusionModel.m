% clear all
% close all

SeanData = readmatrix('/home/omeacock/Documents/SPRruns/patchVelocityRuns/hit_counts_Density_kill1000long_corrected.csv');

numSims = 20;
tSteps = 201;
dt = 5;
w = 1; %Cell width

%Values for the diffusion model
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell 
dx = 10; %Spatial resolution of the coarse-grained grid
noX = 20;
noY = 20;
xWidth = noX * dx;
yHeight = noY * dx;
noReps = 1; %Number of times you should repeat the simulation with different starting conditions to get your confidence interval

noRho = 4;
noLam = 1;
noForce = 5;
noAtFrac = 1;

tList = linspace(0,(tSteps-1)*dt,tSteps);

timecourses = SeanData(1:tSteps*numSims,8)/3699;
lams = SeanData(1:tSteps:tSteps*numSims,5);
forces = SeanData(1:tSteps:tSteps*numSims,2);
rho0s = SeanData(1:tSteps:tSteps*numSims,3);
atFracs = SeanData(1:tSteps:tSteps*numSims,4)/100;
vels = SeanData(1:tSteps:tSteps*numSims,9);

lams = reshape(lams,noAtFrac,noRho,noForce,noLam);
forces = reshape(forces,noAtFrac,noRho,noForce,noLam);
rho0s = reshape(rho0s,noAtFrac,noRho,noForce,noLam);
atFracs = zeros(noAtFrac,noRho,noForce,noLam); %Will load from actual data (input parameter doesn't exactly match the simulation output)
timecourses = reshape(timecourses,tSteps,noAtFrac,noRho,noForce,noLam);
vels = reshape(vels,noAtFrac,noRho,noForce,noLam);

[forces,sortInds] = sort(forces,3);
sortInds = squeeze(sortInds(1,1,:,1));
lams = lams(:,:,sortInds,:);
rho0s = rho0s(:,:,sortInds,:);
timecourses = timecourses(:,:,:,sortInds,:);
vels = vels(:,:,sortInds,:);

figure

for atFracInd = 1:noAtFrac
    for lamInd = 1:noLam
        for fInd = 1:noForce
            for rhoInd = 1:noRho
                subplot(2,2,noRho-rhoInd+1)
                hold on
                ax = gca;
                
                l = lams(atFracInd,rhoInd,fInd,lamInd);
                v = vels(atFracInd,rhoInd,fInd,lamInd);
                r = rho0s(atFracInd,rhoInd,fInd,lamInd);
                a = 1-timecourses(1,atFracInd,rhoInd,fInd,lamInd);
                atFracs(atFracInd,rhoInd,fInd,lamInd) = a;

                lineCol = [(fInd-1)/(noForce-1),1-(fInd-1)/(noForce-1),1];
                
                traceStore = zeros(tSteps,noReps);
                
                for rep = 1:noReps
                    %Initial conditions
                    [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,r,a);
                    pops = cat(3,startA,startS,zeros(noX,noY,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

                    %Allow the system to equilibrate contact compartments without killing or
                    %diffusion
                    [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
                    pops = pops(end,:);
                    pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
                    popsTcourse = pops;

                    %Outer loop - macroscopic mixing (diffusion)
                    for t = 1:(tSteps-1)
                        %Run diffusion of each of the populations separately (can we do this?
                        %Will they be guaranteed to add up to one everywhere if they're run
                        %separately?)
                        for i = 1:size(pops,3)
                            pops(:,:,i) = diffTimestepCN(pops(:,:,i),dt,dx,alphaD*v,true);
                        end

                        %Inner loop - microscopic mixing (contact swapping)
                        [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,l,noX,noY,noConts,noHitBins),[0,dt],pops(:));
                        pops = pops(end,:);
                        pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
                        popsTcourse = cat(4,popsTcourse,pops);
                    end
                    traceStore(:,rep) = squeeze(sum(sum(sum(popsTcourse(:,:,2:6:37,:),3),2),1))/sum(startS(:) + startA(:));
                end
                
                if fInd == 1
                    lH(rhoInd) = plot(tList,squeeze(timecourses(:,atFracInd,rhoInd,fInd,lamInd)),'Color',lineCol);
                else
                    plot(tList,squeeze(timecourses(:,atFracInd,rhoInd,fInd,lamInd)),'Color',lineCol)
                end
                
                plotStdAsArea(tList,traceStore,lineCol,ax)

                pause(0.1)
            end
        end
    end
end

for rhoInd = 1:noRho
    subplot(2,2,rhoInd)
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    axis(ax,[0,1000,0,1])
    title(ax,['\rho_0 = ',num2str(rho0s(1,noRho-rhoInd+1,1,1))])
end

%subplot(1,4,1)
%legend(lH,{'\rho_0 = 0.00025','\rho_0 = 0.001','\rho_0 = 0.004','\rho_0 = 0.016'},'Location','SouthWest')