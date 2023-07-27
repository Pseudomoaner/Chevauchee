function [homoRate,contExRate] = calcMixingContributions(v,lam,noX,noY,noConts,noHitBins,diffDt,dx,D,pops)

%Initial, unperturbed hit rate
hitRate0 = calcHitRate(pops,noHitBins,noConts,noX,noY,lam);

%Apply homogenization

%Run diffusion of each of the populations separately
popsHomo = zeros(size(pops));
for i = 1:size(pops,3)
    popsHomo(:,:,i) = diffTimestepCN(pops(:,:,i),diffDt,dx,D,true);
end

%Need to equilibrate the contact number distributions
popsHomo = equilibrateContNos(popsHomo,noX,noY,noHitBins,noConts);

%Calculate homogenized hit rate
hitRateHomo = calcHitRate(popsHomo,noHitBins,noConts,noX,noY,lam);

%Apply contact exchange
[~,popsContEx] = ode45(@(t,y)diffusiveODEs(t,y,v,0,noX,noY,noConts,noHitBins),[0,diffDt],popsHomo(:));
popsContEx = reshape(popsContEx(end,:),noY,noX,noHitBins*(noConts+1) + 1);

%Calculate homogenization plus contact exchange hit rate
hitRateContEx = calcHitRate(popsContEx,noHitBins,noConts,noX,noY,lam);

%Calculate the hit rate contributions
homoRate = (hitRateHomo - hitRate0)/diffDt;
contExRate = (hitRateContEx - hitRateHomo)/diffDt;
end


% Function that finds the steady-state contact number distribution shape
% for a given attacker fraction and redistributes populations to match
function redistPops = equilibrateContNos(pops,noX,noY,noHitBins,noConts)

redistPops = zeros(size(pops));

for i = 1:noY
    for j = 1:noX
        popsBin = reshape(pops(i,j,2:end),[noHitBins,noConts+1]);

        currDist = sum(popsBin,1);
        ssDist = binopdf(0:noConts,noConts,pops(i,j,1))*(1-pops(i,j,1)); %Steady-state distribution of contact nos
        distDiff = ssDist-currDist;
        
        lostPopFracs = popsBin(:,distDiff<0)./repmat(sum(popsBin(:,distDiff<0)),noHitBins,1);
        lostPops = lostPopFracs.*repmat(distDiff(distDiff<0),noHitBins,1);
        gainedPopSizes = -sum(lostPops,2); %The total size of each hit bin to be redistributed among the gaining contact bins
        gainBinFracs = distDiff(distDiff>0)/sum(distDiff(distDiff>0));
        gainedPops = gainedPopSizes.*repmat(gainBinFracs,noHitBins,1);

        %Subtract the lost populations
        popsBin(:,distDiff<0) = popsBin(:,distDiff<0) + lostPops;
        popsBin(:,distDiff>0) = popsBin(:,distDiff>0) + gainedPops;
        
        redistPops(i,j,1) = pops(i,j,1);
        redistPops(i,j,2:end) = popsBin(:);
    end
end

end

% Function that calculates the current average hit rate for the input
% population distribution
function hitRate = calcHitRate(pops,noHitBins,noConts,noX,noY,lam)

hitRates = zeros(noY,noX);

for i = 1:noY
    for j = 1:noX
        popsBin = reshape(pops(i,j,2:end),[noHitBins,noConts+1]);
    
        %Calculate the baseline kill rate
        unhitPops = squeeze(popsBin(1,:));
        hitRates(i,j) = sum(unhitPops.*(0:noConts))*lam;
    end
end

hitRate = mean(hitRates(:));
end