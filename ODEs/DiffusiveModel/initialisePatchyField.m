function [outA,outS] = initialisePatchyField(dx,xWidth,yHeight,seedDensity,seedFrac)
%INITIALISEPATCHYFIELD creates a new continnum lattice with populations of
%attackers and sensitives distributed according to a Voronoi tesselation
%over seed cells at the specified densities.
%
%   INPUTS:
%       - dx: Spacing of the lattice sites
%       - xWidth, yHeight: Size of the lattice grid, in physical units
%       - seedDensity: The number of seed cells per unit area to be applied
%       - seedFrac: The fraction of seed cells that should be members of
%       the attacker population
%   
%   OUTPUTS:
%       - [outA,outS]: The fraction of attacker (outA) and sensitive (outS)
%       cells at each site in the lattice.
%
%   Author: Oliver J. Meacock, 2023

us = 10; %Upsampling factor for high-resolution grid
tol = 0.05; %Fractional tolerance for the final proportion of attackers to sensitives (must be within x percent of seedFrac)

%Create fine-resolution mesh that will be assigned to seeds and then
%coarse-grained to get final resolution verion of attacker and sensitive
%grids.
[xMesh,yMesh] = meshgrid(dx/(us*2):dx/us:xWidth-dx/(us*2),dx/(us*2):dx/us:yHeight-dx/(us*2));

noSeeds = round(xWidth*yHeight*seedDensity);
noSeeds1 = round(noSeeds * seedFrac);
noSeeds2 = noSeeds - noSeeds1;
seedsX = rand(noSeeds,1)*xWidth;
seedsY = rand(noSeeds,1)*yHeight;
tmpSeeds = [zeros(noSeeds2,1);ones(noSeeds1,1)];
popLabels = tmpSeeds(randperm(noSeeds)); %0 indicates pop 2, 1 pop 1

if noSeeds1 == 0 && seedFrac ~= 0
    error('Seeding density is too low for any attackers to be initialized. Try increasing domain size.')
end

%Create a tiling of these initial seeds so the resulting Voronoi
%diagram is periodic
seedsX = [seedsX;seedsX;seedsX;seedsX - xWidth;seedsX - xWidth;seedsX - xWidth;seedsX + xWidth;seedsX + xWidth;seedsX + xWidth];
seedsY = [seedsY;seedsY - yHeight;seedsY + yHeight;seedsY;seedsY - yHeight;seedsY + yHeight;seedsY;seedsY - yHeight;seedsY + yHeight];
popLabels = repmat(popLabels,9,1);

DT = delaunayTriangulation([seedsX,seedsY]);
nearestSeeds = nearestNeighbor(DT,xMesh(:),yMesh(:));

noPop1 = sum(popLabels(nearestSeeds));
realFrac = noPop1/numel(nearestSeeds);

tooMany = realFrac > seedFrac + (seedFrac * tol);
tooFew = realFrac < seedFrac - (seedFrac * tol);

while tooMany || tooFew
    seedsX = rand(noSeeds,1)*xWidth;
    seedsY = rand(noSeeds,1)*yHeight;

    popLabels = tmpSeeds(randperm(noSeeds)); %0 indicates pop 2, 1 pop 1

    seedsX = [seedsX;seedsX;seedsX;seedsX - xWidth;seedsX - xWidth;seedsX - xWidth;seedsX + xWidth;seedsX + xWidth;seedsX + xWidth];
    seedsY = [seedsY;seedsY - yHeight;seedsY + yHeight;seedsY;seedsY - yHeight;seedsY + yHeight;seedsY;seedsY - yHeight;seedsY + yHeight];
    popLabels = repmat(popLabels,9,1);

    DT = delaunayTriangulation([seedsX,seedsY]);
    nearestSeeds = nearestNeighbor(DT,xMesh(:),yMesh(:));

    noPop1 = sum(popLabels(nearestSeeds));
    realFrac = noPop1/numel(nearestSeeds);

    tooMany = realFrac > seedFrac + (seedFrac * tol);
    tooFew = realFrac < seedFrac - (seedFrac * tol);
end

%Assign population labels to mesh and coarse-grain back to original
%resolution
popA = reshape(popLabels(nearestSeeds) == 1,size(xMesh,1),size(xMesh,2));
popS = reshape(popLabels(nearestSeeds) == 0,size(xMesh,1),size(xMesh,2));

outA = zeros(size(popA,1)/us,size(popA,2)/us);
outS = zeros(size(popA,1)/us,size(popA,2)/us);

for i = 1:size(popA,1)/us
    for j = 1:size(popA,2)/us
        outA(i,j) = mean(mean(popA((i-1)*us + 1:i*us,(j-1)*us + 1:j*us)));
        outS(i,j) = mean(mean(popS((i-1)*us + 1:i*us,(j-1)*us + 1:j*us)));
    end
end