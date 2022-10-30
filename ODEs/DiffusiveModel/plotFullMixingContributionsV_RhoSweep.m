clear all
close all

lam = 0.01;
atFrac = 1/10;

vs = [1,0.1,0.01];
rho0s = [1,0.1,0.01,0.001];

dx = 10; %Granularity of coarse-grained lattice
xWidth = 200;
yHeight = 200;
noX = xWidth/dx;
noY = yHeight/dx;

tMax = 1000;
diffDt = 5; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt + 1;

noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

figure

for vInd = 1:size(vs,2)
    v = vs(vInd);

    baseStore = zeros(size(rho0s));
    contExStore = zeros(size(rho0s));
    homoStore = zeros(size(rho0s));

    for rhoInd = 1:size(rho0s,2)
        rho0 = rho0s(rhoInd);
        
        vList = v*ones(noDiffTsteps,1);
        [baseStore(rhoInd),contExStore(rhoInd),homoStore(rhoInd)] = calcFullMixingContributions(dx,xWidth,yHeight,rho0,atFrac,noConts,noHitBins,noX,noY,lam,vList,diffDt,tMax);
    end
    subplot(size(vs,2),1,vInd)
    bar([baseStore;contExStore;homoStore]','stacked')
    hold on
    plot([0,numel(rho0s)+1],1-[atFrac,atFrac],'k--','LineWidth',1.5)
    xlabel('Seeding density')
    ylabel('Contribution to killing')
    axis([0.4,4.6,0,1])
    title(['v = ',num2str(v)])
    ax = gca;
    ax.XTickLabel = cellstr(num2str(rho0s(:)));
end

legend('Base','Contact exchange','Homogenization')