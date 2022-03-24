function simTracks = simulateCellTracks(fieldsFlow,x,y,seedX,seedY,dt,wrapSets)
%SIMULATEPASSIVETRACKS simulates the trajectories of passive particles in a
%flowing system, given the flowfields.
%
%   INPUTS:
%       -fieldsFlow: The flowfields of your system, either extracted from
%       coarsegrained agent-based simulations or directly from continuum
%       simulations. Should be on a regular lattice. NxMxTx2 array.
%       -x,y: The coordinates of the lattice points. NxM arrays
%       -startX,startY: The coordinates of the seeding points. Column vectors.
%       -wrapSets: Settings for defining wrap-around (if applicable). Set
%       wrapSets.active to true to apply periodic boundary conditions, then
%       use wrapSets.maxX, wrapSets.maxY to specify the boundaries of the
%       periodic domain.
%
%   OUTPUTS:
%       -simTracks: The tracks of the simulated passive objects, formatted
%       using the FAST standard. Must contain at least the x, y, phi, 
%       theta, vmag and times fields.
%
%   Author: Oliver J. Meacock, 2021

%Step through set number of timepoints in the simulated fields, and interpolate
%the velocity and orientation fields. Will use nearest-neighbour for the
%time being for the orientation fields to avoid the circular wrap around
%issues

%Preallocate tracks
for i = 1:size(seedX,1)*noTgts
    simTracks(i).x = zeros(trackLen,1);
    simTracks(i).y = zeros(trackLen,1);
    simTracks(i).theta = zeros(trackLen-1,1);
    simTracks(i).vmag = zeros(trackLen-1,1);
end

%Periodise the x and y grids if necessary, as well as the flowfields
if wrapSets.active
    x = repmat([x-wrapSets.maxX,x,x+wrapSets.maxX],3,1);
    y = repmat([y-wrapSets.maxY;y;y+wrapSets.maxY],1,3);
    fieldsFlow = repmat(fieldsFlow,3,3,1,1);
end

for s = 1:size(seedX,1)
    currX = seedX(s);
    currY = seedY(s);
    
    for i = 1:size(currX,1)
        simTracks(s).x(1) = currX(i);
        simTracks(s).y(1) = currY(i);
        simTracks(s).times = 1:maxT;
    end
    
    for t = 1:maxT
        currFlowU = fieldsFlow(:,:,t,1);
        currFlowV = fieldsFlow(:,:,t,2);

        uQ = interp2(x,y,currFlowU,currX,currY);
        vQ = interp2(x,y,currFlowV,currX,currY);
        
        dX = uQ*dt;
        dY = vQ*dt;
        vmags = sqrt(uQ.^2 + vQ.^2);
        thets = atan2d(dY,dX);
        currX = currX + dX;
        currY = currY + dY;

        %Apply periodic correction to updated positions, if necessary
        if wrapSets.active
            currX = mod(currX,wrapSets.maxX);
            currY = mod(currY,wrapSets.maxY);
        end
        
        for i = 1:size(currX,1)
            simTracks(i+offset).x(t+1) = currX(i);
            simTracks(i+offset).y(t+1) = currY(i);
            simTracks(i+offset).theta(t) = thets(i);
            simTracks(i+offset).vmag(t) = vmags(i);
        end
    end
end

%Invert measurements to bring into line with standard for target tracks.
for i = 1:size(simTracks,2)
    simTracks(i).phi = -simTracks(i).phi;
    simTracks(i).theta = -simTracks(i).theta;
end