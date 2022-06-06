function dydt = wellMixedODEs(t,y,v,lam,noConts,noHitBins)

alphE = 0.11946; %Encounter rate proportionality constant
re = v*alphE;

dydt = zeros(size(y));

J = y(1);
y = y(2:end);

for c = 0:noConts %sensitive-attacker contact index
    cInd = c*noHitBins;
    for h = 1:noHitBins %hit index
        if c == 0
            encounterTerm = re * ((c+1) * y(cInd + h + noHitBins) * (1-J) - (noConts - c) * y(cInd + h) * J);
        elseif c == noConts
            encounterTerm = re * ((noConts + 1 - c) * y(cInd + h - noHitBins) * J - c * y(cInd + h) * (1-J));
        else
            encounterTerm = re * ((c+1) * y(cInd + h + noHitBins) * (1-J) + (noConts + 1 - c) * y(cInd + h - noHitBins) * J - y(cInd + h) * (c * (1-J) + (noConts-c) * J));
        end

        if h == 1
            hitTerm = c * lam * (- y(cInd + h));
        elseif h == noHitBins
            hitTerm = c * lam * (y(cInd + h - 1));
        else
            hitTerm = c * lam * (y(cInd + h - 1) - y(cInd + h));
        end

        dydt(cInd + h + 1) = encounterTerm + hitTerm;
    end
end