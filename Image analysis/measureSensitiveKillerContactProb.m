function pK = measureSensitiveKillerContactProb(imgS,imgK,pxSize)

% l = 2.6; %Average cell length (um)
% w = 1.0; %Average cell width (um)
% 
% aByL = ((l-w)*w + pi*(w/2)^2)/(pi*w + (l-w)*2); %Area to perimeter ratio for an average cell (units = um)

aByL = 0.4313; %Area to perimeter ratio for an average cell (units = um). Empirically derived - see AreaPerimeterMeasurements.xlsx for measurements

boundImg = bwmorph(and(imdilate(imgS,strel('disk',2)),imdilate(imgK,strel('disk',2))),'skel',Inf);

L = sum(boundImg(:)) * pxSize;
As = sum(imgS(:)) * (pxSize^2);

pK = aByL * (L/As);