function [GFPseg,RFPseg] = splitFluo(BFseg,GFPraw,RFPraw,GFPflat,RFPflat)

se = strel('disk',2); %Used for performing the morphological closure on the initial segmentation (gets rid of feathery bits)

%% Option 1: Gaussian mixture model
% GFPsmooth = imgaussfilt(GFPraw./GFPflat,5);
% RFPsmooth = imgaussfilt(RFPraw./RFPflat,5);
% 
% fluoRatSmooth = log(RFPsmooth./GFPsmooth);
% 
% data = fluoRatSmooth(logical(BFseg));
% if numel(data) > 1000000 %This prevents fitgmdist from running out of memory for very large images
%     data = data(randperm(numel(data),1000000));
% end
% model = fitgmdist(data,2); %Mixed Gaussian model of ratiometric pixel distribution in first frame
% idx = cluster(model,data);
% cluster1 = data(idx == 1);
% cluster2 = data(idx == 2);
% fluoThresh = min(max(cluster1),max(cluster2));

%% Option 2: split by mean flatfield-corrected background
% fluoThresh = 0;

%% Option 3: assume 50:50 population distribution
GFPsmooth = imgaussfilt(GFPraw./GFPflat,5);
RFPsmooth = imgaussfilt(RFPraw./RFPflat,5);

fluoRatSmooth = log(RFPsmooth./GFPsmooth);

data = fluoRatSmooth(logical(BFseg));

fluoThresh = prctile(data,50);

%Once you have the fluorescence threshold from the smoothed data,
%return to less smoothed data for higher resolution results
fluoRat = log(imgaussfilt(RFPraw./RFPflat,2)./imgaussfilt(GFPraw./GFPflat,2));

GFPsegTmp = and(fluoRat <= fluoThresh, BFseg);
RFPsegTmp = and(fluoRat > fluoThresh, BFseg);

GFPseg = imopen(GFPsegTmp,se);
RFPseg = imopen(RFPsegTmp,se);