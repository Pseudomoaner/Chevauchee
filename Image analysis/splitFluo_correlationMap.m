function [GFPseg,RFPseg] = splitFluo(CMseg,GFPraw,RFPraw,GFPflat,RFPflat)

se = strel('disk',2); %Used for performing the morphological closure on the initial segmentation (gets rid of feathery bits)

%In the first part, we will use the fluorescence data to refine the
%brightfield segmentation (if the field isn't already confluent)
if mean(CMseg(:)) < 0.8 %I.e. if non-confluent
    fluoSum = RFPraw./RFPflat + GFPraw./GFPflat;
    imgRc = (fluoSum-prctile(fluoSum(:),0.1))./(prctile(fluoSum(:),99.9)-prctile(fluoSum(:),0.1)); %Rescaled image between 0 and 1
    imgRc = imgaussfilt(imgRc,2);
    T = adaptthresh(imgRc,0.7);
    CMseg = and(imbinarize(imgRc,T),CMseg);
end

%Now calculate the threshold for splitting your two populations
GFPsmooth = imgaussfilt(GFPraw./GFPflat,5);
RFPsmooth = imgaussfilt(RFPraw./RFPflat,5);

fluoRatSmooth = log(RFPsmooth./GFPsmooth);

%% Option 1: Gaussian mixture model
if sum(CMseg(:)) > 1
    data = fluoRatSmooth(logical(CMseg));
    if numel(data) > 1000000 %This prevents fitgmdist from running out of memory for very large images
        data = data(randperm(numel(data),1000000));
    end
    model = fitgmdist(data,2); %Mixed Gaussian model of ratiometric pixel distribution in first frame

    %If either model doesn't converge or component fractions are way off,
    %revert to default threshold of 50th percentile
    if ~model.Converged || sum(model.ComponentProportion < 0.4) > 0
        fluoThresh = prctile(data,50);
    else
        idx = cluster(model,data);
        cluster1 = data(idx == 1);
        cluster2 = data(idx == 2);
        fluoThresh = min(max(cluster1),max(cluster2));
    end
else %If there is no material detected, both GFP and RFP segmentations will be empty - so choose an arbitrary threshold to prevent crashing
    fluoThresh = 0;
end

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

%Once you have the fluorescence threshold from the smoothed data,
%return to less smoothed data for higher resolution results
fluoRat = log(imgaussfilt(RFPraw./RFPflat,2)./imgaussfilt(GFPraw./GFPflat,2));

GFPsegTmp = and(fluoRat <= fluoThresh, CMseg);
RFPsegTmp = and(fluoRat > fluoThresh, CMseg);

GFPseg = imopen(GFPsegTmp,se);
RFPseg = imopen(RFPsegTmp,se);
