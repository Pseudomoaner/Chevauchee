clear all
close all

root = '/home/omeacock/Documents/SPRruns/Driveby/HitTypeChecks';
branch = 'SimulationResults_F=%f_HT=%s.mat';

forces = [0.5,0.75,1,1.25];

noConts = 5; %Average number of contacts per cell

figure
hold on
ax = gca;

for fInd = 1:size(forces,2)
    f=forces(fInd);
    
    load(fullfile(root,sprintf(branch,f,'constant')))
    
    tSeries = 0:fieldSettings.dt:fieldSettings.dt*fieldSettings.maxF;
    
    noHit = cellfun(@(x)sum(x>0),trackableData.Hit);
    noAt = sum(trackableData.FireRate{1} > 0);
    noSens = numel(trackableData.FireRate{1}) - noAt;
        
    constTrace = (noSens-noHit)/(noSens+noAt);
    
    load(fullfile(root,sprintf(branch,f,'distributed')))
    
    noHit = cellfun(@(x)sum(x>0),trackableData.Hit);
    noAt = sum(trackableData.FireRate{1} > 0);
    noSens = numel(trackableData.FireRate{1}) - noAt;
        
    distTrace = (noSens-noHit)/(noSens+noAt);
    
    lineCol = [(fInd-1)/(size(forces,2)-1),1-(fInd-1)/(size(forces,2)-1),1];
        
    plot(ax,tSeries,constTrace,'Color',lineCol)
    plot(ax,tSeries/noConts,distTrace,'--','Color',lineCol)
end
    