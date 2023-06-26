clear all
close all

lam = 0.01; %Firing rate
atFrac = 1/10;

vs = [0.03,0.1,0.3];
rho0s = [0.256,0.064,0.016,0.004,0.001,0.00025];

dx = 10; %Granularity of coarse-grained lattice
xWidth = 300;
yHeight = 300;
noX = xWidth/dx;
noY = yHeight/dx;

tMax = 1000;
diffDt = 5; %Timestep between diffusion timesteps
noDiffTsteps = tMax/diffDt + 1;

noHitBins = 2; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

figure
cList = {[198, 159, 137]/255,[206, 45, 45]/255,[119, 225, 119]/255,[245,245,245]/255};

for rhoInd = 1:size(rho0s,2)
    rho0 = rho0s(rhoInd);

    baseStore = zeros(size(vs));
    contExStore = zeros(size(vs));
    homoStore = zeros(size(vs));

    [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);

    for vInd = 1:size(vs,2)
        v = vs(vInd);
        
        vList = v*ones(noDiffTsteps,1);
        [baseStore(vInd),contExStore(vInd),homoStore(vInd)] = calcFullMixingContributions(dx,startA,startS,atFrac,noConts,noHitBins,noX,noY,lam,vList,diffDt,tMax,1);
        
        subplot(size(vs,2),size(rho0s,2),(vInd-1)*size(rho0s,2)+rhoInd)
        donut([baseStore(vInd),contExStore(vInd),homoStore(vInd),1-atFrac-baseStore(vInd)-contExStore(vInd)-homoStore(vInd)],[],cList)
    end
%     subplot(1,size(rho0s,2),rhoInd)
%     br = bar([baseStore;contExStore;homoStore]','stacked','LineWidth',1);
%     hold on
%     for k = 1:3
%         br(k).FaceColor = cList(k,:);
%     end
%     plot([0,numel(rho0s)+1],1-[atFrac,atFrac],'k--','LineWidth',1.5)
%     xlabel('Average cell velocity')
%     ylabel('Contribution to killing')
%     axis([0.4,3.6,0,1])
%     title(['\rho_0 = ',num2str(rho0)])
%     ax = gca;
%     ax.XTickLabel = cellstr(num2str(vs(:)));
%     ax.LineWidth = 1.5;
end

% legend('Base','Contact exchange','Homogenization')