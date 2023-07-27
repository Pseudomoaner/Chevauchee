function [outA,outS] = initialisePatchyFieldSpecified(dx,xWidth,yHeight,patchSpecs)
%INITIALISEPATCHYFIELDSPECIFIED is used to specify a field of Voronoi-type
%patches if the specifications for the seeding points are already known.
%Typically this is because they have already been generated for the SPR
%model.
%
%   INPUTS:
%       - dx: Spacing of the lattice sites
%       - xWidth, yHeight: Size of the lattice grid, in physical units
%       - patchSpecs: A structure containing the locations and identities
%       of the seed cells.
%   
%   OUTPUTS:
%       - [outA,outS]: The fraction of attacker (outA) and sensitive (outS)
%       cells at each site in the lattice.
%
%   Author: Oliver J. Meacock, 2023

seedsX = patchSpecs.seedsX;
seedsY = patchSpecs.seedsY;
popLabels = patchSpecs.popLabels;

us = 10; %Upsampling factor for high-resolution grid

%Create fine-resolution mesh that will be assigned to seeds and then
%coarse-grained to get final resolution verion of attacker and sensitive
%grids.
[xMesh,yMesh] = meshgrid(dx/(us*2):dx/us:xWidth-dx/(us*2),dx/(us*2):dx/us:yHeight-dx/(us*2));

DT = delaunayTriangulation([seedsX,seedsY]);
nearestSeeds = nearestNeighbor(DT,xMesh(:),yMesh(:));

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