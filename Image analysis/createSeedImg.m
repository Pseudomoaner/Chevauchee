clear all
close all

root = 'C:\Users\omeacock\Desktop\SeanData\';
branches = {'OD1','OD0_1','OD0_01','OD0_001'};
atTwig = '_AttackSeedsImg.tif';
sensTwig = '_SensSeedsImg.tif';
outTwig = '_Voroni.tif';

atC = [77,153,229]/255;
sensC = [255,127,12]/255;

dx = 1;

for i = 1:size(branches,2)
    atImg = logical(imread([root,branches{i},atTwig]));
    sensImg = logical(imread([root,branches{i},sensTwig]));

    %Weed out any tiny regions that made it through manual correction
    atSeedImg = bwpropfilt(atImg,'Area',[5,inf]);
    sensSeedImg = bwpropfilt(sensImg,'Area',[5,inf]);

    atSeeds = regionprops(atSeedImg,'Centroid');
    sensSeeds = regionprops(sensSeedImg,'Centroid');
    atSeedPos = vertcat(atSeeds.Centroid);
    sensSeedPos = vertcat(sensSeeds.Centroid);

    xWidth = size(atSeedImg,2);
    yHeight = size(atSeedImg,1);

    [xMesh,yMesh] = meshgrid(dx/2:dx:xWidth-dx/2,dx/2:dx:yHeight-dx/2);

    seedPos = [atSeedPos;sensSeedPos];
    popLabels = [zeros(size(atSeedPos,1),1);ones(size(sensSeedPos,1),1)];

    DT = delaunayTriangulation([seedPos(:,1),seedPos(:,2)]);
    nearestSeeds = nearestNeighbor(DT,xMesh(:),yMesh(:));

    popA = reshape(popLabels(nearestSeeds) == 0,size(xMesh,1),size(xMesh,2));
    popS = reshape(popLabels(nearestSeeds) == 1,size(xMesh,1),size(xMesh,2));
    
    rCh = ((popA + atSeedImg) * atC(1))/2 + ((popS + sensSeedImg) * sensC(1))/2;
    gCh = ((popA + atSeedImg) * atC(2))/2 + ((popS + sensSeedImg) * sensC(2))/2;
    bCh = ((popA + atSeedImg) * atC(3))/2 + ((popS + sensSeedImg) * sensC(3))/2;

    outImg = cat(3,rCh,gCh,bCh);

    imwrite(outImg,[root,branches{i},outTwig])
end