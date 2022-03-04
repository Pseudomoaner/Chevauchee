function dydt = patchyODEs(t,y,v,lam,rho0,tgtJ)

alphE = 0.641; %Encounter rate proportionality constant
alphM = 1.0163; %Mixing rate proportionality constant

noBins = (numel(y)-1)/2;
dydt = zeros(size(y));

J = y(end);

%Do the changing 'attacked' bin size
dydt(end) = sqrt(rho0)*alphM*v*(tgtJ - J);

%Do the 'protected' bins
for i = 1:noBins
    dydt(i) = v*alphE*((1-J)*y(i+noBins) - J*y(i)) - dydt(end) * y(i);
end

%Do the 'attacked' bins
for i = 1:noBins
    if i == 1
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) - lam*(y(i+noBins)) + dydt(end) * y(i);
    elseif i == noBins
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) + lam*(y(i+noBins-1)) + dydt(end) * y(i);
    else
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) + lam*(y(i+noBins-1) - y(i+noBins)) + dydt(end) * y(i);
    end
end