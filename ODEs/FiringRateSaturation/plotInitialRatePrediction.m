clear all
close all

root = '/home/omeacock/Documents/SPRruns/VelocitySaturationRuns';
branch = 'SimulationResults_lam_%f_force_%f.mat';

forces = [0,0.5,0.75,1,1.25,1.5,1.75,2];
lams = [0.005,0.01,0.02,0.05];

vs = 0:0.001:0.3;
atFrac = 0.01;

alphE = 0.1033; %Per-contact exchange rate proportionality constant
noConts = 5; %Average number of contacts per cell

figure
hold on
ax = gca;

vMeans = zeros(size(forces,2),size(lams,2));
kFRs = zeros(size(forces,2),size(lams,2));

%Plot the simulation datapoints
for lamInd = 1:size(lams,2)
    lam = lams(lamInd);
    for fInd = 1:size(forces,2)
        f = forces(fInd);
        
        load(fullfile(root,sprintf(branch,lam,f)))
        
        initConts = sum(endField.fireCells > 0) * noConts; %The number of hit cells you can discount because they formed part of the initial transient      
        vMeans(fInd,lamInd) = mean(arrayfun(@(x)mean(x.vmag),data));
        
        %Ignore times 1-200, as these will represent an ititial transient
        %as attackers hit their immediate neighbours. So use only times
        %200-500
        hitNoChange = sum(trackableData.Hit{end}>0) - initConts;
        includeTime = fieldSettings.maxF * fieldSettings.dt;
        killRate = hitNoChange/includeTime; %Number of sensitives newly hit per unit time
        kFRs(fInd,lamInd) = killRate/numel(endField.xCells); %Fraction of sensitives newly hit per unit time
    end
end

%Plot the theory lines
for lamInd = 1:size(lams,2)
    lam = lams(lamInd);

    re = alphE * vs;
    predRates = lam * atFrac * (re./(re*(2-atFrac)+lam)) * noConts;
 
    plot(vs,predRates,'Color',[(lamInd-1)/3,1-(lamInd-1)/3,(lamInd-1)/3],'LineWidth',1.5)
    plot(vMeans(:,lamInd),kFRs(:,lamInd),'o','MarkerEdgeColor',[(lamInd-1)/3,1-(lamInd-1)/3,(lamInd-1)/3],'MarkerFaceColor',[0.8,0.8,0.8],'LineWidth',1.5)
end

xlabel('Velocity')
ylabel('Initial killing rate')

ax.Box = 'on';
ax.LineWidth = 1.5;
