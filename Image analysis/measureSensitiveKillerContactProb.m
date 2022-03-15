function pK = measureSensitiveKillerContactProb(imgS,imgK,pxSize)

aByL = 0.245; %Area to perimeter ratio for an average cell (units = um)

boundImg = bwmorph(and(imdilate(imgS,strel('disk',2)),imgK),'skel',Inf);

L = sum(boundImg(:)) * pxSize;
As = sum(imgS(:)) * (pxSize^2);

pK = aByL * (L/As);