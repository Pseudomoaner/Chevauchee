clear all
close all

rho0s = [1,0.1,0.01,0.001]*0.064;
colours = [10,58,92;31,126,193;76,178,250;156,212,252]/255;

lams = [0.001,0.002,0.005,0.01,0.02]; %CDI firing rate (s^-1) - 0.001 corresponds to 3.6 firings / hr
hitEfficiencies = [0.01,0.02,0.05,0.1,0.2]; %The impact of each hit on the targeted cell's ongoing growth rate
atFrac = 5.5/11; %Attacker fraction
noReps = 5;

CIstore = zeros(size(rho0s,2),noReps);

for lamInd = 1:size(lams,2)
    lam = lams(lamInd);
    for efInd = 1:size(hitEfficiencies,2)
        hitEfficiency = hitEfficiencies(efInd);
        
        subplot(size(lams,2),size(hitEfficiencies,2),(efInd-1)*size(lams,2)+lamInd)
        ax = gca;
        hold(ax,'on')
        
        for i = 1:4
            rho0 = rho0s(i); %Seeding density (cells cellWidth^-2)
            for r = 1:noReps
                load(sprintf('C:\\Users\\Olivier\\OneDrive - Université de Lausanne\\Simulations\\Driveby\\popDat_rho_%f_ef_%f_lam_%f_r_%i.mat',rho0,hitEfficiency,lam,r))
            
                CIstore(i,r) = CIs(end);
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
        axis(ax,[0.0005,5,1,1e27])
        xticks(1:4)
        yticks([1,1e10,1e20])
        title(['\lambda = ', num2str(lam), ', e = ', num2str(hitEfficiency)])
    end
end
