function dydt = patchyODEsVariableV(t,y,vmax,vrate,confT,lam,rho0,tgtJ)
%vmax = maximum velocity
%vrate = width of velocity Gaussian
%confT = confluency time

alphE = 0.083; %Encounter rate proportionality constant
alphM = 0.79472; %Mixing rate proportionality constant

v = vmax*exp(-((t-confT)/vrate)^2);

re = v*alphE;

noBins = (numel(y)-1)/5;
dydt = zeros(size(y));

J = y(end);

%Do the changing 'attacked' bin size
dydt(end) = sqrt(rho0)*alphM*v*(tgtJ - J);

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