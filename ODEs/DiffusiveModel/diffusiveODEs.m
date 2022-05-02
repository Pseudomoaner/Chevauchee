function dydt = diffusiveODEs(t,y,v,lam,noX,noY,noConts,noHitBins)

alphE = 0.11946; %Encounter rate proportionality constant
re = v*alphE;

%Reshape y to recreate the attacker and sensitive arrays
popBlock = reshape(y,noY,noX,noHitBins*(noConts+1) + 1);
popBlockSens = popBlock(:,:,2:end);

%J is equal to the sensitive-attacker contact probability - the fraction of
%attackers within a given bin
J = popBlock(:,:,1);

updatePopBlock = zeros(size(popBlockSens));

for c = 0:noConts %sensitive-attacker contact index
    cInd = c*noHitBins;
    for h = 1:noHitBins %hit index
        if c == 0
            encounterTerm = re * ((c+1) * popBlockSens(:,:,cInd + h + noHitBins) .* (1-J) - (noConts-c) * popBlockSens(:,:,cInd + h) .* J);
        elseif c == noConts
            encounterTerm = re * ((noConts + 1 -c) * popBlockSens(:,:,cInd + h - noHitBins) .* J - c * popBlockSens(:,:,cInd + h) .* (1-J));
        else
            encounterTerm = re * ((c+1) * popBlockSens(:,:,cInd + h + noHitBins) .* (1-J) + (noConts + 1 -c) * popBlockSens(:,:,cInd + h - noHitBins) .* J - popBlockSens(:,:,cInd + h) .* (c * (1-J) + (noConts-c) * J));
        end

        if h == 1
            hitTerm = c * lam * (-popBlockSens(:,:,cInd + h));
        elseif h == noHitBins
            hitTerm = c * lam * (popBlockSens(:,:,cInd + h - 1));
        else
            hitTerm = c * lam * (popBlockSens(:,:,cInd + h - 1) - popBlockSens(:,:,cInd + h));
        end

        updatePopBlock(:,:,cInd + h) = encounterTerm + hitTerm;
    end
end

dydt = cat(3,zeros(noY,noX),updatePopBlock);
dydt = dydt(:);