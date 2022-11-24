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
noDiffTsteps = tMax/diffDt;
alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate

noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
noConts = 5; %Number of contacts made by each cell

figure

for vInd = 1:size(vs,2)
    v = vs(vInd);
    D = v*alphaD;
    for rhoInd = 1:size(rho0s,2)
        rho0 = rho0s(rhoInd);

        %Initial conditions
        [startA,startS] = initialisePatchyField(dx,xWidth,yHeight,rho0,atFrac);
        pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

        %Allow the system to equilibrate contact compartments without killing or
        %diffusion
        [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
        pops = pops(end,:);
        pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
        popsTcourse = pops;

        %Storage for mixing type contributions
        homoStore = zeros(noDiffTsteps,1);
        contExStore = zeros(noDiffTsteps,1);

        %Outer loop - macroscopic mixing (diffusion)
        for t = 1:noDiffTsteps
            %Calculate the relative contributions of homogenization and contact
            %exchange, respectively
            [homoStore(t),contExStore(t)] = calcMixingContributions(v,lam,noX,noY,noConts,noHitBins,diffDt,dx,D,pops);

            %Run diffusion of each of the populations separately
            for i = 1:size(pops,3)
                pops(:,:,i) = diffTimestepCN(pops(:,:,i),diffDt,dx,D,true);
            end

            %Inner loop - microscopic mixing (contact swapping)
            [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,lam,noX,noY,noConts,noHitBins),[0,diffDt],pops(:));
            pops = pops(end,:);
            pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
            popsTcourse = cat(4,popsTcourse,pops);
        end

        subplot(size(vs,2),size(rho0s,2),(vInd-1)*size(rho0s,2) + rhoInd)
        plot(0:diffDt:(tMax-diffDt),homoStore,'r','LineWidth',1.5)
        hold on
        plot(0:diffDt:(tMax-diffDt),contExStore,'b','LineWidth',1.5)
        ax = gca;
        ax.LineWidth = 1.5;
        %axis(ax,[0,1000,0,3e-4])
        title(ax,['\rho_0 = ',num2str(rho0),', v = ',num2str(v)])
    end
end