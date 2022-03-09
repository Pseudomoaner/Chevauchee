function dydt = wellMixedODEs(t,y,v,lam,J)

alphE = 0.083; %Encounter rate proportionality constant
re = v*alphE;

noBins = numel(y)/5;
dydt = zeros(size(y));

for c = 0:4 %sensitive-attacker contact index
    cInd = c*noBins;
    for h = 1:noBins %hit index
        if c == 0
            encounterTerm = re * ((c+1) * y(cInd + h + noBins) * (1-J) - (4-c) * y(cInd + h) * J);
        elseif c == 4
            encounterTerm = re * ((5-c) * y(cInd + h - noBins) * J - c * y(cInd + h) * (1-J));
        else
            encounterTerm = re * ((c+1) * y(cInd + h + noBins) * (1-J) + (5-c) * y(cInd + h - noBins) * J - y(cInd + h) * (c * (1-J) + (4-c) * J));
        end

        if h == 1
            hitTerm = c * lam * (- y(cInd + h));
        elseif h == noBins
            hitTerm = c * lam * (y(cInd + h - 1));
        else
            hitTerm = c * lam * (y(cInd + h - 1) - y(cInd + h));
        end

        dydt(cInd + h) = encounterTerm + hitTerm;
    end
end