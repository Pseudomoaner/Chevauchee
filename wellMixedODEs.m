function dydt = wellMixedODEs(t,y,v,lam)

alphE = 0.641; %Encounter rate proportionality constant

noBins = numel(y)/2;
dydt = zeros(size(y));

J = sum(y(noBins+1:end))/sum(y);

%Do the 'protected' bins
for i = 1:noBins
    dydt(i) = v*alphE*((1-J)*y(i+noBins) - J*y(i));
end

%Do the 'attacked' bins
for i = 1:noBins
    if i == 1
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) - lam*(y(i+noBins));
    elseif i == noBins
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) + lam*(y(i+noBins-1));
    else
        dydt(i+noBins) = v*(J*y(i) - (1-J)*y(i+noBins)) + lam*(y(i+noBins-1) - y(i+noBins));
    end
end