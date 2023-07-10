%Code to generate IBM and continuum simulations for Fig. 1.
%Note you need both QuasiS3R and move-n-murder on your path - unlike other
%QuasiS3R runs, this one also does appropriately-matched continuum
%simulations.
clear all
close all

%% File control and parameter definitions

%Define file names etc.
RootSim = '/home/omeacock/Documents/SPRruns/patchVelocityRuns';
outputSprMatName = 'SprResults_f_%f_rho0_%f.mat';
outputContinuumMatName = 'ContinuumResults_f_%f_rho0_%f.mat';

%Settings applied to entire field
fieldSettings.xWidth = 200; %Width of the simulated domain (in units of fieldSettings.lam)
fieldSettings.yHeight = 200; %Height of the simulated domain (in units of fieldSettings.lam)
fieldSettings.postDivMovement = 'reverse'; %How the daughter cell should move following cell division. Either 'reverse' (the opposite direction to the mother) or 'same' (the same direction as the mother).
fieldSettings.fireRange = 2; %Range of the CDI system
fieldSettings.killType = 'husk'; %Whether to remove cells from simulation after death ('lyse') or inactivate them, but leave their bodies ('husk')
fieldSettings.killThresh = 1000; %Number of hits needed to kill a cell
fieldSettings.hitRateType = 'constant'; %Whether CDI hits are diluted over all cells ('distributed') or the per-neighbour hit rate is kept the same regardless of the number of contacts ('constant')
fieldSettings.growthRate = 0.0; %Average increase in aspect ratio over one unit of time.
fieldSettings.divThresh = 8; %Aspect ratio at which the cell should divide.
fieldSettings.areaFrac = 0.35; %Fraction of the total area that should be occupied by cells.

%Choose your cell and barrier settings
barrierSettings.type = 'none'; %Type of static barriers that should be present in simulation - either none or loaded

%Settings for the active rods - note that the use of the 'LatticedXYCells'
%option means that all rods are assumed to be identical. Initialization can
%be customized by writing additional code in the WensinkField.populateField
%function.
cellSettings.type = 'LatticedXYCells'; %Type of rod initialization conditions that should be applied - either singleCell, doubleCell or LatticedXYCells
cellSettings.a = 4; %Aspect ratio of rods (relative to fieldSettings.lam)
cellSettings.r = 0; %Reversal rate associated with each rod
cellSettings.c = [1,0,0]; %RGB values for the colour you want to make the cells of population 1
cellSettings.f = 1;
cellSettings.fire = 0; %Firing rate of CDI system for this population
cellSettings.pop = 's'; %Population label to specify which cells can kill each other.

%Settings for the patch initialization, which runs after initial active
%configuration has been reached
patchSettings.patchType = 'Voronoi';
patchSettings.seedFrac = 0.1; % Fraction of seeds that should be associated with population 2
patchSettings.reversalRate = 0;
patchSettings.colour = [0,0.5,1];
patchSettings.fireRate = 0.02;
patchSettings.popLabel = 't';
patchSettings.force = 1;
patchSettings.tol = 0.05; %Tolerance for actual population fraction (real value guaranteed between seedFrac or fraction - for 'Voronoi' and 'Circle' patch types respectively - plus/minus fraction * tol)

%Output settings
dispSettings.saveFrames = false; %Whether or not to save visualisations of each sampled timepoint
dispSettings.ImgPath = 'Frame_%04d.tif'; %Generic name for each output frame (will be fed into sprintf, so use appropriate string formatting)
dispSettings.colourCells = 'Hits'; %How rods should be recoloured at each sampling point. If set to 'None', will retain any previously set colour.
dispSettings.saveType = 'draw'; %Type of method used to visualise rods - either 'plot' or 'draw'. 'plot' will produce and save a Matlab figure, while 'draw' will draw ellipses directly into an image.
dispSettings.imagedirectory = [RootSim,filesep,'ColourCells']; %Defines where the output images will be located
if ~exist(dispSettings.imagedirectory,'dir') %Set up visualisation directory
    mkdir(dispSettings.imagedirectory);
end

%Processing settings
procSettings.startTime = 1; %Can be useful to cut out some times, if it takes some time for the model to settle down
procSettings.velocitySmoothingSize = 2; %Size of the smoothing window used to smooth rod position data
procSettings.minTrackLength = 5; %Minimum length of a track to be kept following track assembly
procSettings.pixSize = 0.2; %In the same units as lam. Value is defined by the settings in WensinkField.drawField() (could bring those parameters out to here in the future)

%Global simulation settings (defined separately from e.g. field settings so
%they can easily applied uniformly during parameter sweeps).
burnInDt = 0.0625;
startMotileDt = 0.1; %Size of the motility timestep (to begin with)
firingDt = 1; %Size of the timestep for calculating firing events
samplingRate = 5.0; %How frequently samples of the simulation should be taken
burnInSimTime = 0;
settlingSimTime = 200; %How long it will take for the simulation to settle into an active configuration
targetSimTime = 1000; %Target motile simulation time
contactFind = false; %Whether or not to return structures containing instantaneous cell-cell contact data

rho0s = [0.256,0.064];
fs = [0,0.5,1,1.5,2];

for rhoInd = 1:size(rho0s,2)
    for fInd = 1:size(fs,2)
        f = fs(fInd);
        rho0 = rho0s(rhoInd);

        patchSettings.seedDensity = rho0;

        %% Part 0: Initialize field for this simulation (including burn-in)
        startField = WensinkField(fieldSettings);
        startField = startField.populateField(barrierSettings,cellSettings,fieldSettings.areaFrac);

        tmpFieldSettings = fieldSettings;
        tmpFieldSettings.burnIndt = burnInDt;
        tmpFieldSettings.burnInSteps = round(burnInSimTime/burnInDt);
        tmpFieldSettings.FrameSkip = round(samplingRate/burnInDt);
        tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
        startField = simulateWensinkFieldBurnIn(startField,tmpFieldSettings,dispSettings);

        %% Part 1: Make sure that the selected value of motiledt doesn't make the simulation explode, or reduce until you reach numerical stability

        %Make sure that this value of motiledt doesn't make the simulation explode.
        tmpFieldSettings = fieldSettings;
        tmpFieldSettings.motiledt = startMotileDt;
        tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
        tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
        [tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
        simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
        simOK = true;
        while ~simOK %Repeatedly reduce simulation step size if needed
            tmpFieldSettings.motiledt = tmpFieldSettings.motiledt;
            tmpFieldSettings.FrameSkip = round(samplingRate/tmpFieldSettings.motiledt);
            tmpFieldSettings.motileSteps = tmpFieldSettings.FrameSkip;
            [tmpPCs,~] = simulateWensinkFieldInitial(startField,tmpFieldSettings,dispSettings);
            simOK = checkTimestepAppropriate(tmpPCs,tmpFieldSettings);
        end
        fieldSettings.motiledt = tmpFieldSettings.motiledt;
        fieldSettings.FrameSkip = round(samplingRate/fieldSettings.motiledt);
        fieldSettings.FireSkip = round(firingDt/fieldSettings.motiledt);

        %% Part 2: Do initial simulation to allow system to reach an active configuration
        fieldSettings.motileSteps = ceil(settlingSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
        [~,intermediateField] = simulateWensinkFieldInitial(startField,fieldSettings,dispSettings);

        %Setup a patch in the centre of the domain containing the second population
        [intermediateField,patchSpecs] = makePatch(intermediateField,patchSettings);
        
        intermediateField.fCells = ones(size(intermediateField.xCells)) * f;

        %% Part 3: Do another (fully sampled) simulation for a longer period of time - only data from this simulation period will be stored
        fieldSettings.motileSteps = ceil(targetSimTime/(fieldSettings.motiledt*fieldSettings.FrameSkip))*fieldSettings.FrameSkip;
        [PCs,endField,hitNos,contactSet] = simulateWensinkField(intermediateField,fieldSettings,dispSettings,contactFind);

        %% Part 4: Process data and save simulation results
        fieldSettings.dt = fieldSettings.motiledt * fieldSettings.FrameSkip;
        fieldSettings.maxF = round(fieldSettings.motileSteps/fieldSettings.FrameSkip);

        [data,trackableData,toMappings,fromMappings] = processModelPCs(PCs,procSettings,fieldSettings);

        fullSprMatOut = fullfile(RootSim,sprintf(outputSprMatName,f,rho0));
        save(fullSprMatOut,'data','trackableData','toMappings','fromMappings','fieldSettings','cellSettings','procSettings','samplingRate','startMotileDt','hitNos','endField','contactSet')
    
        %% Part 5: Run equivalent continuum model with same input parameters and starting configuration
        dX = 10; %Size of the coarse-grained grid (in cell widths)
        alphaD = 5.638; %Proportionality constant that converts velocity into cell diffusion rate
        noHitBins = 6; %Number of different hit categories to take into consideration, from 0 to noHitBins-1 plus
        noConts = 5; %Average number of contacts made by each cell to neighbours
        tSteps = targetSimTime/fieldSettings.dt;

        v = mean(arrayfun(@(x)mean(x.vmag),data));

        noX = round(fieldSettings.xWidth/dX);
        noY = round(fieldSettings.yHeight/dX);

        [startA,startS] = initialisePatchyFieldSpecified(dX,fieldSettings.xWidth,fieldSettings.yHeight,patchSpecs);
        pops = cat(3,startA,startS,zeros(noY,noX,noHitBins*(noConts+1)-1)); %Create population array, attackers in first layer, unhit sensitives in second)

        %Allow the system to equilibrate contact compartments without killing or
        %diffusion
        [t,pops] = ode45(@(t,y)diffusiveODEs(t,y,1,0,noX,noY,noConts,noHitBins),[0,100],pops(:));
        pops = pops(end,:);
        pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
        popsTcourse = pops;

        %Outer loop - macroscopic mixing (diffusion)
        for t = 1:tSteps
            %Run diffusion of each of the populations separately
            for i = 1:size(pops,3)
                pops(:,:,i) = diffTimestepCN(pops(:,:,i),fieldSettings.dt,dX,alphaD*v,true);
            end

            %Inner loop - microscopic mixing (contact swapping)
            [~,pops] = ode45(@(t,y)diffusiveODEs(t,y,v,patchSettings.fireRate,noX,noY,noConts,noHitBins),[0,fieldSettings.dt],pops(:));
            pops = pops(end,:);
            pops = reshape(pops,noY,noX,noHitBins*(noConts+1) + 1);
            popsTcourse = cat(4,popsTcourse,pops);
        end

        fullContinuumMatOut = fullfile(RootSim,sprintf(outputContinuumMatName,f,rho0));
        save(fullContinuumMatOut, 'popsTcourse', 'patchSpecs')
    end
end
