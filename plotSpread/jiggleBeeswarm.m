function jiggledX = jiggleBeeswarm(rawX,rawY,binWidth)
%JIGGLEBEESWARM creates a copy of the input dataset with non-overlapping
%datapoints. Based on the plotSpread.m fuction of https://uk.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot?s_tid=mwa_osa_a

stdWidth = 1;

if isempty(binWidth) %Implies you want to spread points based on their relative (rather than absolute) closeness in X
    binWidth = 0.3;
    
    %% TRANSFORM DATA
    % Here, I try to estimate what the aspect ratio of the data is going to be
    nData = size(rawX,1);
    fh = figure('Visible','off');
    minMax = [min(rawY);max(rawY)];
    
    plot([0.5;nData+0.5],minMax,'o');
    
    aspectRatio = get(gca,'DataAspectRatio');
    close(fh);
    
    tFact = aspectRatio(2)/aspectRatio(1);
    
    %% SPREAD POINTS
    % assign either nData, or xValues number of values, in case we're working
    % with group-indices
    
    % transform and sort
    currentData = rawY / tFact;
    
    % add x
    currentData = [ones(size(rawY))*rawX,currentData];
else
    currentData = [ones(size(rawY))*rawX,rawY];
end


% step through the data in 0.1 increments. If there are multiple
% entries, spread along x
for y = min(currentData(:,2)):binWidth:max(currentData(:,2))
    % find values
    valIdx = find(currentData(:,2) >= y & currentData(:,2) < y+binWidth);
    
    nVal = length(valIdx);
    if nVal > 1
        % spread
        spreadWidth = stdWidth*0.9*(1-exp(log(0.9)*(nVal-1)));
        spreadDist = spreadWidth / (nVal - 1);
        if rem(nVal,2) == 0
            offset = spreadDist / 2;
        else
            offset = eps;
        end
        for v = 1:nVal
            currentData(valIdx(v),1) = rawX + offset;
            % update offset
            offset = offset - sign(offset) * spreadDist * v;
        end
    end
end

jiggledX = currentData(:,1);
