function mixVar = measureSensitiveKillerIntermixing(GFPseg,RFPseg,pxSize)

boxEdgeSize = 10*0.8; %Size of coarse-grained grid in um (here 10 cell widths multiplied by the average 0.8 um cell width)
boxEdgeSizePx = boxEdgeSize/pxSize; %Size of coarse-grained grid in pixels

minVol = 0.05; %Minimum fraction a given grid element must be filled to be included in calculation

noBoxX = round(size(GFPseg,2)/boxEdgeSizePx);
noBoxY = round(size(GFPseg,1)/boxEdgeSizePx);

xList = round(linspace(1,size(GFPseg,2),noBoxX+1));
yList = round(linspace(1,size(GFPseg,1),noBoxY+1));

coarseG = zeros(noBoxY,noBoxX);
coarseR = zeros(noBoxY,noBoxX);
boxA = zeros(noBoxY,noBoxX);

for i = 1:noBoxX
    for j = 1:noBoxY
        Gchunk = GFPseg(yList(j):yList(j+1),xList(i):xList(i+1));
        Rchunk = RFPseg(yList(j):yList(j+1),xList(i):xList(i+1)); 

        coarseG(j,i) = sum(Gchunk(:));
        coarseR(j,i) = sum(Rchunk(:));
        boxA(j,i) = numel(Gchunk);
    end
end

totFill = coarseR + coarseG; %Total number of pixels filled by either population
badLocs = (totFill./boxA) < minVol;

mixVar = var(coarseR(~badLocs)./(coarseR(~badLocs)+coarseG(~badLocs)));