function fluoThresh = findFluoThresh(BFseg,GFPraw,RFPraw,GFPflat,RFPflat)

se = strel('disk',2); %Used for performing the morphological closure on the initial segmentation (gets rid of feathery bits)

%In the first part, we will use the fluorescence data to refine the
%brightfield segmentation (if the field isn't already confluent)
if mean(BFseg(:)) < 0.8 %I.e. if non-confluent
    fluoSum = RFPraw./RFPflat + GFPraw./GFPflat;
    imgRc = (fluoSum-prctile(fluoSum(:),0.1))./(prctile(fluoSum(:),99.9)-prctile(fluoSum(:),0.1)); %Rescaled image between 0 and 1
    imgRc = imgaussfilt(imgRc,2);
    T = adaptthresh(imgRc,0.7);
    BFseg = and(imbinarize(imgRc,T),BFseg);
end

%Now calculate the threshold for splitting your two populations
GFPsmooth = imgaussfilt(GFPraw./GFPflat,2);
RFPsmooth = imgaussfilt(RFPraw./RFPflat,2);

fluoRatSmooth = log(RFPsmooth./GFPsmooth);

%% Option 1: Gaussian mixture model
% if sum(BFseg(:)) > 1
%     data = fluoRatSmooth(logical(BFseg));
%     if numel(data) > 1000000 %This prevents fitgmdist from running out of memory for very large images
%         data = data(randperm(numel(data),1000000));
%     end
%     model = fitgmdist(data,2); %Mixed Gaussian model of ratiometric pixel distribution in first frame
% 
%     %If either model doesn't converge or component fractions are way off,
%     %revert to default threshold of 50th percentile
%     if ~model.Converged || sum(model.ComponentProportion < 0.4) > 0
%         fluoThresh = prctile(data,50);
%     else
%         idx = cluster(model,data);
%         cluster1 = data(idx == 1);
%         cluster2 = data(idx == 2);
%         fluoThresh = min(max(cluster1),max(cluster2));
%     end
% else %If there is no material detected, both GFP and RFP segmentations will be empty - so choose an arbitrary threshold to prevent crashing
%     fluoThresh = 0;
% end

%% Option 2: split by mean flatfield-corrected background
% fluoThresh = 0;

%% Option 3: assume 50:50 population distribution
% GFPsmooth = imgaussfilt(GFPraw./GFPflat,5);
% RFPsmooth = imgaussfilt(RFPraw./RFPflat,5);
% 
% fluoRatSmooth = log(RFPsmooth./GFPsmooth);
% 
% data = fluoRatSmooth(logical(BFseg));
% 
% fluoThresh = prctile(data,50);

%% Option 4: split by 0-gradient point in pdf closest to 50:50 in the cdf
bw = 0.01;
fracThresh = 0.2;

if sum(BFseg(:)) > 1
    data = fluoRatSmooth(logical(BFseg));
    [N,E] = histcounts(data,'BinWidth',bw);
    
    %Chop off the tails of the distribution
    Ncdf = cumsum(N)/sum(N);
    N = N(and(Ncdf < 0.99,Ncdf > 0.01));
    cents = E([and(Ncdf < 0.99,Ncdf > 0.01),false]) + bw/2; %Bin centers
    
    minima = find(diff(sign(diff(N)))==2) + 1; %Locations of local minima
    
    %If a minimum can't be found, assume 50:50 distribution
    %as backup.
    if isempty(minima)
        fluoThresh = prctile(data,50);
    elseif numel(minima) == 1
        fluoThresh = cents(minima);
    elseif numel(minima) == 2 %If you have two minima, choose the average position between as your threshold
        fluoThresh = cents(round(mean(minima)));
    else  %If you have many minima, select the one with the lowest associated number of pixels
        [~,fluoThreshLoc] = min(N(minima));
        fluoThresh = cents(minima(fluoThreshLoc));
    end
    
    fracA = sum(N(cents < fluoThresh))/sum(N);
    
    %If an intermediate minimum can't be found, assume 50:50 distribution
    %as backup.
    if fracA < fracThresh || fracA > (1-fracThresh)
        fluoThresh = prctile(data,50);
    end
else %If there is no material detected, both GFP and RFP segmentations will be empty - so choose an arbitrary threshold to prevent crashing
    fluoThresh = 0;
end