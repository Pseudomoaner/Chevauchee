clear all
close all

lams = [0.005,0.01,0.02,0.05];
vs = 0:0.01:1;
atFrac = 0.01;

alphE = 0.1;

figure
hold on
ax = gca;

for lamInd = 1:size(lams,2)
    lam = lams(lamInd);

    re = alphE * vs;
    predRates = lam * atFrac * (re./(re*(2-atFrac)+lam));

    plot(vs,predRates)
end

xlabel('Velocity')
ylabel('Initial killing rate')