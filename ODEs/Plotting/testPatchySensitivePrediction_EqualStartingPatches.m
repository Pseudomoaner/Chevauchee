%Plots the results of the QuasiS3R runs produced by
%runQuasiS3rsweep_PropulsionForce_Patchy.m along with continuum simulations
%initialised with the same seed distribution.

clear all
close all

Root = '/home/omeacock/Documents/SPRruns/patchVelocityRuns';
SPRbranch = 'SprResults_f_%f_rho0_%f.mat';
ContBranch = 'ContinuumResults_f_%f_rho0_%f.mat';

rhos = [0.256,0.064,0.016,0.004,0.001,0.00025];
fs = 0:0.5:2;

dX = 10; %Size of the coarse-grained grid (in cell widths)
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Average number of contacts made by each cell to neighbours
fireRate = 0.02;

figure

for rhoInd = 1:size(rhos,2)
    subplot(1,size(rhos,2),rhoInd)
    ax = gca;
    hold(ax,'on');
    
    rho = rhos(rhoInd);
    for fInd = 1:size(fs,2)
        f = fs(fInd);
        
        load(fullfile(Root,sprintf(SPRbranch,f,rho)))
        load(fullfile(Root,sprintf(ContBranch,f,rho)))
        
        tSeries = 0:fieldSettings.dt:fieldSettings.dt*fieldSettings.maxF;
        
        v = mean(arrayfun(@(x)mean(x.vmag),data));
        
%         noX = round(fieldSettings.xWidth/dX);
%         noY = round(fieldSettings.yHeight/dX);
%         
%         [startA,startS] = initialisePatchyFieldSpecified(dX,fieldSettings.xWidth,fieldSettings.yHeight,patchSpecs);
%         pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)
%         
%         %Allow the system to equilibrate contact compartments without killing or
%         %diffusion
%         [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
%         pops = pops(end,:);
%         pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
%         popsTcourse = pops;
%         
%         %Outer loop - macroscopic mixing (diffusion)
%         for t = 1:size(tSeries,2)-1
%             %Run diffusion of each of the populations separately
%             for i = 1:size(pops,3)
%                 pops(:,:,i) = diffTimestepCN(pops(:,:,i),fieldSettings.dt,dX,alphaD*v,true);
%             end
% 
%             %Inner loop - microscopic mixing (contact swapping)
%             [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,fireRate,noX,noY,noConts,noHitBins),[0,fieldSettings.dt],pops(:));
%             pops = pops(end,:);
%             pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
%             popsTcourse = cat(4,popsTcourse,pops);
%         end
        
        contTrace = squeeze(sum(sum(sum(popsTcourse(:,:,2:6:37,:),3),2),1))/(size(popsTcourse,1)*size(popsTcourse,2));
        
        noHit = cellfun(@(x)sum(x>0),trackableData.Hit);
        noAt = sum(trackableData.FireRate{1} > 0);
        noSens = numel(trackableData.FireRate{1}) - noAt;
        
        sprTrace = (noSens-noHit)/(noSens+noAt);
        
        lineCol = [(fInd-1)/(size(fs,2)-1),1-(fInd-1)/(size(fs,2)-1),1];
        
        plot(ax,tSeries,sprTrace,'Color',lineCol)
        plot(ax,tSeries,contTrace,'--','Color',lineCol)
    end
end

for rhoInd = 1:size(rhos,2)
    subplot(1,size(rhos,2),rhoInd)
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    axis(ax,[0,1000,0,1])
    title(ax,['\rho_0 = ',num2str(rhos(rhoInd))])
end

%subplot(1,4,1)
%legend(lH,{'\rho_0 = 0.00025','\rho_0 = 0.001','\rho_0 = 0.004','\rho_0 = 0.016'},'Location','SouthWest')